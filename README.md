# GAME-CpG-Predictor

We designed the CpG Predictor to serve as a baseline for models, as CpG content is correlated with many functional readouts. It ignores the type, cell type, and species keys in the Evaluator requests since it is simply counting CpG instances across the sequences. For point requests, it will return the mean CpG counts across the input sequence, and, for track requests, it will return the CpG density at each base pair (calculated using a sliding window of 50bp). Any scaling type can be fulfilled by this Predictor.

## Important Links

- To learn more about the GAME Framework: [Main GAME Repository](https://github.com/de-Boer-Lab/Genomic-API-for-Model-Evaluation), [preprint](https://www.biorxiv.org/content/10.1101/2025.07.04.663250v1.full)
- GAME Documentation: [ReadTheDocs](https://genomic-api-for-model-evaluation-documentation.readthedocs.io)
- Pre-built CpG Predictor container image: [Hugging Face](https://huggingface.co/datasets/deBoerLab/CpG_Predictor_GAME)
- List of all [GAME Modules](https://github.com/de-Boer-Lab/GAME_modules)

## Features

- **Point Predictions**: Returns average CpG count across entire sequence
- **Track Predictions**: Returns per-base CpG density using sliding window approach
- **Flexible Scaling**: Linear or log2-scaled outputs (default is linear)
- **Multiple Formats**: Supports JSON and MessagePack request/response formats
- **Containerized**: Apptainer container for reproducible deployment

### Run CpG Prediction

```bash
apptainer run --containall cpg_predictor.sif <ip_address> <port>
# Example: apptainer run --containall cpg_predictor.sif 0.0.0.0 5000
```

### Build Container

If you would like to make edits to the code/re-build the container use the command below. Run the build from inside `src/` (the `.def` copies files from the current directory).

```bash
cd src
apptainer build cpg_predictor.sif cpg_predictor.def
```

## API Endpoints

### GET `/help`

Returns predictor metadata and version information.

### GET `/formats`

Returns supported request/response formats (JSON, MessagePack).

### POST `/predict`

Main prediction endpoint. Accepts DNA sequences and returns CpG predictions.

**Request Format:**

```json
{
  "readout": "point" # or "track",
  "prediction_tasks": [
    {
      "name": "task_name",
      "type": "experiment_type",
      "cell_type": "cell_type",
      "species": "species",
      "scale": "linear" # or "log"
    }
  ],
  "sequences": {
    "seq_id": "ACGTACGT..."
  },
  "prediction_ranges": {
    "seq_id": [start, end]  # 0-based, half-open [start, end): start < end <= len(sequence)
  }
}
```

**Response Format:**

```json
{
  "predictor_name": "CpG_20261029-234547_EDT",  # CpG_<build timestamp>
  "bin_size": 1, # only returned for track requests
  "prediction_tasks": [
    {
      "name": "task_name",
      "type_requested": "accessibility",
      "type_actual": ["cpg_density"],
      "cell_type_requested": "cell_type",
      "cell_type_actual": null,  # cell type is ignored by this predictor
      "species_requested": "mus_musculus",
      "species_actual": null,  # species is ignored by this predictor
      "scale_prediction_requested": "linear",
      "scale_prediction_actual": "linear",
      "predictions": {
        "seq_id": [0.025, 0.032, ...]
      }
    }
  ]
}
```

**Errors:**

- **400** `bad_prediction_request`: missing mandatory keys, invalid values, or malformed `prediction_ranges` (not a list, negative indices, or `start >= end`)
- **422** `prediction_request_failed`: a `readout` other than `point`/`track`, or a `prediction_ranges` entry out of bounds for its sequence

## Algorithm Details

- **Point Mode**: Calculates mean CpG frequency across entire sequence, as CpGs per bp: `(CpG count + 1e-9) / sequence length`. Log scale returns `log2` of this value (the `1e-9` avoids `log2(0)`)
- **Track Mode**: Uses 50bp sliding window to compute per-base CpG density
- **CpG Definition**: Dinucleotide sequence "CG"
- **Density Calculation**: (CpG count / window length) × 100
- **Prediction ranges** (half-open `[start, end)`):
  - Point: the CpG mean is calculated over `sequence[start:end]` only
  - Track: returns `end - start` values, identical to computing the density over the full sequence and keeping positions `start` to `end - 1` (windows near the range edges still see the flanking sequence)

## Directory Structure

```
src/
├── cpg_predictor.def            # Apptainer build definition
├── config.py                    # Predictor name, formats, paths
├── predictor_RestAPI.py         # Main Flask API server
├── cpg_utils.py                 # CpG calculation functions
├── predictor_help_message.json  # Metadata
├── error_checking_functions.py  # Validation utilities
├── schema_validation.py         # Request schema validation
└── predictor_content_handler.py # Format encoding/decoding
```
