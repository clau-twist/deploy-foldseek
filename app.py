import os
import subprocess
import pandas as pd
import logging
import shutil
import tempfile

from flask import Flask, jsonify, request
from typing import Dict, Any

# Set up basic logging
logger = logging.getLogger(__name__)

app = Flask(__name__)

@app.route('/list-files')
def list_files():
    directory="/home/ec2-user/databases"
    try:
        files = os.listdir(directory)
        # Filter out only files (not subdirectories)
        files = [f for f in files if os.path.isfile(os.path.join(directory, f))]
        return jsonify(files)
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route("/foldseek", methods=["post"])
def run_foldseek(): # The arguments are removed from the function signature
    """
    Performs the FoldSeek search.
    It writes the input PDB to a temporary file and executes the FoldSeek binary.
    """
    if not request.is_json:
        return jsonify({"error": "Request must be in JSON format"}), 400

    available_databases = ["afdb-sp", "afdb-up", "esmdb"]
    data = request.get_json()
    query_pdb = data.get("query_pdb")
    database = data.get("database")

    if not query_pdb or database is None:
        return jsonify({"error": "Missing 'query_pdb' or 'database' in request body"}), 400
    if database not in available_databases:
        return jsonify({"error": f"Database is not an available database. Choose from {available_databases}"})

    logger.info("Starting FoldSeek prediction...")

    temp_dir = tempfile.mkdtemp()
    try:
        with tempfile.NamedTemporaryFile(mode="w+", suffix=".pdb", delete=True) as temp_pdb_file:
            # Define paths for temporary files
            query_pdb_path = os.path.join(temp_dir, "query.pdb")
            results_path = os.path.join(temp_dir, "results.tsv")
            tmp_work_path = os.path.join(temp_dir, "tmp")

            with open(query_pdb_path, "w") as f:
                f.write(query_pdb)

            headers = [
                "mismatch", "gapopen", "qstart", "qend", "tstart", "tend",
                "evalue", "prob", "qlen", "tlen", "qaln", "taln",
                "qseq", "tseq"
            ]
            if database != "esmdb":
                headers.append("taxname")
            command = [
                "foldseek",
                "easy-search",
                query_pdb_path,
                f"/home/ec2-user/databases/{database}",
                results_path,
                tmp_work_path,
                "--format-mode", "4",
                "--format-output",
                "mismatch,gapopen,qstart,qend,tstart,tend,evalue,prob,qlen,tlen,qaln,taln,qseq,tseq,taxname",
            ]

            logger.info(f"Running command: {' '.join(command)}")

            process = subprocess.run(
                command,
                capture_output=True,
                text=True,
                check=True  # This will raise an exception on non-zero exit codes
            )

            logger.info(f"FoldSeek stdout: {process.stdout}")

            df = pd.read_csv(
                results_path,
                sep='\t',
                names=headers,
                header=None,
                skip_blank_lines=True
            )
            csv_results = df.to_csv(index=False)

            logger.info("FoldSeek search completed successfully.")
            return jsonify({csv_results}), 200

    except subprocess.CalledProcessError as e:
        logger.error(f"FoldSeek execution failed with exit code {e.returncode}")
        logger.error(f"FoldSeek stderr: {e.stderr}")
        error_code = jsonify({f"Foldseek stderr: {e.stderr}"}), 400
    finally:
        shutil.rmtree(temp_dir)
    return error_code

def get_top_n_results(tsv_path: str, top_n: int = 10) -> Dict[str, Any]:
    """
    Parse TSV results from Foldseek, sort by lowest evalue, and return top N results.
    Uses pandas for efficient parsing and sorting.

    Args:
        tsv_content: TSV content string from Foldseek results
        top_n: Number of top results to return (default: 10)

    Returns:
        Dictionary containing top N results sorted by evalue (lowest first)
    """
    # Column headers for Foldseek TSV output
    headers = [
        "mismatch", "gapopen", "qstart", "qend", "tstart", "tend",
        "evalue", "prob", "qlen", "tlen", "qaln", "taln",
        "qseq", "tseq", "taxname"
    ]

    try:
        # Read TSV data into pandas DataFrame
        df = pd.read_csv(
            tsv_path,
            sep='\t',
            names=headers,
            header=None,
            skip_blank_lines=True
        )

        # Drop rows with missing evalue
        df = df.dropna(subset=['evalue'])

        # Convert numeric columns to appropriate types
        numeric_cols = ["mismatch", "gapopen", "qstart", "qend", "tstart", "tend", "qlen", "tlen"]
        for col in numeric_cols:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors='coerce')

        # Ensure evalue and prob are float
        df['evalue'] = pd.to_numeric(df['evalue'], errors='coerce')
        df['prob'] = pd.to_numeric(df['prob'], errors='coerce')

        # Get top N results with lowest evalue
        df_sorted = df.nsmallest(top_n, 'evalue')

        # Convert to list of dictionaries
        top_results = df_sorted.to_dict('records')

        return {
            "total_matches": len(df),
            "top_n": top_n,
            "results": top_results
        }

    except Exception as e:
        return {"error": f"Error parsing results: {str(e)}", "results": []}


if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000)

