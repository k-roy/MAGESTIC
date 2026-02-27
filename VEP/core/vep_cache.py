#!/usr/bin/env python3
"""
VEP Cache Module - Sequence-Hash-Based Prediction Caching

Stores VEP predictions indexed by sequence hash + model parameters to avoid
redundant computation across jobs and projects.

Usage:
    from vep_cache import VEPCache

    cache = VEPCache(CACHE_DIR, "shorkie")

    # Define model parameters (affects cache key)
    MODEL_PARAMS = {"window_size": 16384, "model_version": "1.0"}
    params_hash = cache.hash_params(MODEL_PARAMS)

    # Check cache before prediction
    seq_hash = cache.hash_sequence(sequence)
    cached = cache.get(seq_hash, params_hash)
    if cached:
        return cached

    # Run prediction and store
    result = model.predict(sequence)
    cache.put(seq_hash, params_hash, result,
              source_project="bloom_qtl", source_test_id=test_id)

Author: Kevin R. Roy
Date: 2026-02-22
"""

import hashlib
import json
import sqlite3
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
from dataclasses import dataclass


@dataclass
class CacheStats:
    """Statistics about cache usage."""
    total_entries: int
    models: List[str]
    oldest_entry: Optional[datetime]
    newest_entry: Optional[datetime]
    total_size_bytes: int


class VEPCache:
    """SQLite-based cache for VEP predictions.

    Indexes predictions by sequence hash + parameter hash, allowing
    reuse across different test IDs that happen to have the same sequence.

    Attributes:
        cache_dir: Base directory for all model caches
        model_name: Name of the model (e.g., "shorkie", "evo2_7b")
        db_path: Path to the SQLite database file
    """

    def __init__(self, cache_dir: Path, model_name: str):
        """Initialize cache for a specific model.

        Args:
            cache_dir: Base directory for VEP cache (e.g., common/vep_cache/)
            model_name: Model name (creates subdirectory if needed)
        """
        self.cache_dir = Path(cache_dir)
        self.model_name = model_name
        self.model_dir = self.cache_dir / model_name
        self.db_path = self.model_dir / "cache.sqlite"

        # Ensure directory exists
        self.model_dir.mkdir(parents=True, exist_ok=True)

        # Initialize database
        self._init_db()

    def _init_db(self) -> None:
        """Initialize the SQLite database with required schema."""
        with sqlite3.connect(str(self.db_path)) as conn:
            cursor = conn.cursor()

            # Enable WAL mode for better concurrent read/write performance
            cursor.execute("PRAGMA journal_mode=WAL")
            cursor.execute("PRAGMA synchronous=NORMAL")
            cursor.execute("PRAGMA cache_size=-64000")  # 64MB cache

            cursor.execute("""
                CREATE TABLE IF NOT EXISTS predictions (
                    seq_hash TEXT NOT NULL,
                    params_hash TEXT NOT NULL,
                    seq_length INTEGER,
                    model_version TEXT,
                    window_size INTEGER,

                    -- Results
                    score REAL,
                    score_json TEXT,

                    -- Metadata
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                    source_project TEXT,
                    source_test_id TEXT,

                    PRIMARY KEY (seq_hash, params_hash)
                )
            """)

            cursor.execute("""
                CREATE INDEX IF NOT EXISTS idx_seq_length
                ON predictions(seq_length)
            """)

            cursor.execute("""
                CREATE INDEX IF NOT EXISTS idx_created
                ON predictions(created_at)
            """)

            cursor.execute("""
                CREATE INDEX IF NOT EXISTS idx_source_project
                ON predictions(source_project)
            """)

            conn.commit()

    @staticmethod
    def hash_sequence(sequence: str) -> str:
        """Generate SHA256 hash of a DNA/protein sequence.

        Args:
            sequence: The sequence to hash

        Returns:
            Hexadecimal hash string (64 characters)
        """
        return hashlib.sha256(sequence.upper().encode()).hexdigest()

    @staticmethod
    def hash_params(params: dict) -> str:
        """Generate MD5 hash of model parameters.

        Args:
            params: Dictionary of model parameters

        Returns:
            Hexadecimal hash string (32 characters)
        """
        # Sort keys for deterministic hashing
        param_str = json.dumps(params, sort_keys=True)
        return hashlib.md5(param_str.encode()).hexdigest()

    def get(self, seq_hash: str, params_hash: str) -> Optional[Dict[str, Any]]:
        """Look up cached prediction.

        Args:
            seq_hash: SHA256 hash of the sequence
            params_hash: MD5 hash of model parameters

        Returns:
            Dictionary with cached results, or None if not found
        """
        with sqlite3.connect(str(self.db_path)) as conn:
            cursor = conn.cursor()

            cursor.execute("""
                SELECT score, score_json, model_version, source_project, source_test_id
                FROM predictions
                WHERE seq_hash = ? AND params_hash = ?
            """, (seq_hash, params_hash))

            row = cursor.fetchone()

        if row is None:
            return None

        score, score_json, model_version, source_project, source_test_id = row

        result = {
            "score": score,
            "model_version": model_version,
            "source_project": source_project,
            "source_test_id": source_test_id,
            "from_cache": True
        }

        # Parse full results if available
        if score_json:
            try:
                result["full_results"] = json.loads(score_json)
            except json.JSONDecodeError:
                pass

        return result

    def put(
        self,
        seq_hash: str,
        params_hash: str,
        result: Dict[str, Any],
        source_project: str,
        source_test_id: str,
        seq_length: Optional[int] = None,
        model_version: Optional[str] = None,
        window_size: Optional[int] = None
    ) -> None:
        """Store prediction in cache.

        Args:
            seq_hash: SHA256 hash of the sequence
            params_hash: MD5 hash of model parameters
            result: Dictionary with prediction results (must have 'score' key)
            source_project: Name of the project that generated this prediction
            source_test_id: Original test_id for reference
            seq_length: Length of the sequence (optional, for filtering)
            model_version: Version of the model used
            window_size: Window size used for prediction
        """
        # Extract primary score
        score = result.get("score") or result.get("logSED") or result.get("delta_logits")

        # Serialize full results
        score_json = json.dumps(result)

        with sqlite3.connect(str(self.db_path)) as conn:
            cursor = conn.cursor()

            cursor.execute("""
                INSERT OR REPLACE INTO predictions
                (seq_hash, params_hash, seq_length, model_version, window_size,
                 score, score_json, source_project, source_test_id)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                seq_hash, params_hash, seq_length, model_version, window_size,
                score, score_json, source_project, source_test_id
            ))

            conn.commit()

    def get_batch(
        self,
        seq_hashes: List[str],
        params_hash: str
    ) -> Dict[str, Optional[Dict[str, Any]]]:
        """Look up multiple cached predictions efficiently.

        Args:
            seq_hashes: List of sequence hashes to look up
            params_hash: MD5 hash of model parameters

        Returns:
            Dictionary mapping seq_hash -> result (or None if not found)
        """
        if not seq_hashes:
            return {}

        with sqlite3.connect(str(self.db_path)) as conn:
            cursor = conn.cursor()

            # Use parameterized query with IN clause
            placeholders = ",".join(["?" for _ in seq_hashes])
            query = f"""
                SELECT seq_hash, score, score_json, model_version, source_project, source_test_id
                FROM predictions
                WHERE params_hash = ? AND seq_hash IN ({placeholders})
            """

            cursor.execute(query, [params_hash] + seq_hashes)
            rows = cursor.fetchall()

        results = {h: None for h in seq_hashes}

        for row in rows:
            seq_hash, score, score_json, model_version, source_project, source_test_id = row
            result = {
                "score": score,
                "model_version": model_version,
                "source_project": source_project,
                "source_test_id": source_test_id,
                "from_cache": True
            }
            if score_json:
                try:
                    result["full_results"] = json.loads(score_json)
                except json.JSONDecodeError:
                    pass
            results[seq_hash] = result

        return results

    def put_batch(
        self,
        entries: List[Tuple[str, Dict[str, Any], str]],
        params_hash: str,
        source_project: str,
        model_version: Optional[str] = None,
        window_size: Optional[int] = None
    ) -> int:
        """Store multiple predictions efficiently.

        Args:
            entries: List of (seq_hash, result_dict, test_id) tuples
            params_hash: MD5 hash of model parameters
            source_project: Name of the project
            model_version: Version of the model used
            window_size: Window size used

        Returns:
            Number of entries stored
        """
        if not entries:
            return 0

        rows = []
        for seq_hash, result, test_id in entries:
            score = result.get("score") or result.get("logSED") or result.get("delta_logits")
            score_json = json.dumps(result)
            seq_length = result.get("seq_length")

            rows.append((
                seq_hash, params_hash, seq_length, model_version, window_size,
                score, score_json, source_project, test_id
            ))

        with sqlite3.connect(str(self.db_path)) as conn:
            cursor = conn.cursor()

            cursor.executemany("""
                INSERT OR REPLACE INTO predictions
                (seq_hash, params_hash, seq_length, model_version, window_size,
                 score, score_json, source_project, source_test_id)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, rows)

            conn.commit()

        return len(rows)

    def get_stats(self) -> CacheStats:
        """Return cache statistics.

        Returns:
            CacheStats with counts and metadata
        """
        with sqlite3.connect(str(self.db_path)) as conn:
            cursor = conn.cursor()

            cursor.execute("SELECT COUNT(*) FROM predictions")
            total_entries = cursor.fetchone()[0]

            cursor.execute("SELECT DISTINCT source_project FROM predictions")
            projects = [row[0] for row in cursor.fetchall() if row[0]]

            cursor.execute("SELECT MIN(created_at), MAX(created_at) FROM predictions")
            oldest, newest = cursor.fetchone()

        # Get database file size
        size_bytes = self.db_path.stat().st_size if self.db_path.exists() else 0

        return CacheStats(
            total_entries=total_entries,
            models=[self.model_name],
            oldest_entry=datetime.fromisoformat(oldest) if oldest else None,
            newest_entry=datetime.fromisoformat(newest) if newest else None,
            total_size_bytes=size_bytes
        )

    def count(self, params_hash: Optional[str] = None) -> int:
        """Count entries in cache.

        Args:
            params_hash: If provided, count only entries with this params_hash

        Returns:
            Number of entries
        """
        with sqlite3.connect(str(self.db_path)) as conn:
            cursor = conn.cursor()

            if params_hash:
                cursor.execute(
                    "SELECT COUNT(*) FROM predictions WHERE params_hash = ?",
                    (params_hash,)
                )
            else:
                cursor.execute("SELECT COUNT(*) FROM predictions")

            count = cursor.fetchone()[0]
        return count

    def clear(self, params_hash: Optional[str] = None) -> int:
        """Clear cache entries.

        Args:
            params_hash: If provided, only clear entries with this params_hash.
                        If None, clears all entries.

        Returns:
            Number of entries deleted
        """
        with sqlite3.connect(str(self.db_path)) as conn:
            cursor = conn.cursor()

            if params_hash:
                cursor.execute(
                    "DELETE FROM predictions WHERE params_hash = ?",
                    (params_hash,)
                )
            else:
                cursor.execute("DELETE FROM predictions")

            deleted = cursor.rowcount
            conn.commit()
        return deleted


def create_model_params_file(
    cache_dir: Path,
    model_name: str,
    params: Dict[str, Any],
    notes: str = ""
) -> None:
    """Create a params.json file documenting model parameters.

    Args:
        cache_dir: Base cache directory
        model_name: Model name
        params: Model parameters dictionary
        notes: Optional notes about the model
    """
    model_dir = cache_dir / model_name
    model_dir.mkdir(parents=True, exist_ok=True)

    params_file = model_dir / "params.json"

    doc = {
        "model_name": model_name,
        "created": datetime.now().isoformat(),
        "parameters": params,
        "notes": notes
    }

    with open(params_file, "w") as f:
        json.dump(doc, f, indent=2)

    print(f"Created {params_file}")


# Convenience function for quick cache lookups
def get_cached_prediction(
    cache_dir: Path,
    model_name: str,
    sequence: str,
    params: Dict[str, Any]
) -> Optional[Dict[str, Any]]:
    """Quick helper to check if a prediction is cached.

    Args:
        cache_dir: Base cache directory
        model_name: Model name (e.g., "shorkie")
        sequence: The DNA/protein sequence
        params: Model parameters

    Returns:
        Cached result or None
    """
    cache = VEPCache(cache_dir, model_name)
    seq_hash = cache.hash_sequence(sequence)
    params_hash = cache.hash_params(params)
    return cache.get(seq_hash, params_hash)


if __name__ == "__main__":
    # Test the cache
    import tempfile

    with tempfile.TemporaryDirectory() as tmpdir:
        cache_dir = Path(tmpdir)
        cache = VEPCache(cache_dir, "test_model")

        # Test hashing
        seq = "ATCGATCGATCG"
        seq_hash = cache.hash_sequence(seq)
        print(f"Sequence hash: {seq_hash[:16]}...")

        params = {"window_size": 16384, "version": "1.0"}
        params_hash = cache.hash_params(params)
        print(f"Params hash: {params_hash}")

        # Test put/get
        result = {"score": 0.5, "logSED": -0.123, "extra": "data"}
        cache.put(seq_hash, params_hash, result,
                  source_project="test", source_test_id="test_1")

        cached = cache.get(seq_hash, params_hash)
        assert cached is not None
        assert cached["score"] == 0.5
        assert cached["from_cache"] is True
        print("Put/get test passed!")

        # Test stats
        stats = cache.get_stats()
        print(f"Cache stats: {stats.total_entries} entries, {stats.total_size_bytes} bytes")

        # Test batch operations
        entries = [
            (cache.hash_sequence(f"SEQ{i}"), {"score": i * 0.1}, f"test_{i}")
            for i in range(10)
        ]
        n_stored = cache.put_batch(entries, params_hash, "batch_test")
        print(f"Stored {n_stored} entries in batch")

        # Test batch get
        hashes = [cache.hash_sequence(f"SEQ{i}") for i in range(10)]
        batch_results = cache.get_batch(hashes, params_hash)
        n_found = sum(1 for v in batch_results.values() if v is not None)
        print(f"Found {n_found}/{len(hashes)} in batch lookup")

        print("\nAll tests passed!")
