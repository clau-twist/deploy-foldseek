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
    # Check if the request contains JSON
    if not request.is_json:
        return jsonify({"error": "Request must be in JSON format"}), 400

    data = request.get_json()
    query_pdb = data.get("query_pdb")
    n = data.get("n")

    # Validate that both parameters were provided
    if not query_pdb or n is None:
        return jsonify({"error": "Missing 'query_pdb' or 'n' in request body"}), 400

    logger.info("Starting FoldSeek prediction...")

    # Create a temporary directory to avoid file collisions between requests
    temp_dir = tempfile.mkdtemp()

    # Slight workaround cuz it's a binary lol
    try:
        with tempfile.NamedTemporaryFile(mode="w+", suffix=".pdb", delete=True) as temp_pdb_file:
            # Define paths for temporary files
            query_pdb_path = os.path.join(temp_dir, "query.pdb")
            results_path = os.path.join(temp_dir, "results.tsv")
            tmp_work_path = os.path.join(temp_dir, "tmp")

            # Write the input PDB string to the temporary file
            with open(query_pdb_path, "w") as f:
                f.write(query_pdb)

            # Construct the command-line arguments for FoldSeek
            # This provides a standard, detailed TSV output
            command = [
                "foldseek",
                "easy-search",
                query_pdb_path,
                "/home/ec2-user/databases/afdb",
                results_path,
                tmp_work_path,
                "--format-mode", "4",
                "--format-output",
                "mismatch,gapopen,qstart,qend,tstart,tend,evalue,prob,qlen,tlen,qaln,taln,qseq,tseq,taxname",
            ]

            logger.info(f"Running command: {' '.join(command)}")

            # Execute the FoldSeek command
            process = subprocess.run(
                command,
                capture_output=True,
                text=True,
                check=True  # This will raise an exception on non-zero exit codes
            )

            logger.info(f"FoldSeek stdout: {process.stdout}")

            results = get_top_n_results(results_path, n)

            # # Read the tab-separated result file
            # with open(results_path, "r") as f:
            #     results = f.read()

            logger.info("FoldSeek search completed successfully.")
            return results

    except subprocess.CalledProcessError as e:
        # Log the specific error from FoldSeek if it fails
        logger.error(f"FoldSeek execution failed with exit code {e.returncode}")
        print(f"FoldSeek execution failed with exit code {e.returncode}")
        logger.error(f"FoldSeek stderr: {e.stderr}")
        print(f"FoldSeek stderr: {e.stderr}")
        raise RuntimeError(f"FoldSeek execution failed: {e.stderr}")
    finally:
        # IMPORTANT: Clean up the temporary directory and its contents
        shutil.rmtree(temp_dir)
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

