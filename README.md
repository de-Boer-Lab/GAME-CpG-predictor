# GAME-CpG-Predictor

We designed the CpG Predictor to serve as a baseline for models, as CpG content is correlated with many functional readouts. It ignores the type, cell type, and species keys in the Evaluator requests since it is simply counting CpG instances across the sequences. For point requests, it will return the mean CpG counts across the input sequence, and, for track requests, it will return the CpG density at each base pair (calculated using a sliding window of 100bp). Any scaling type can be fulfilled by this Predictor. 

## Features

- **Point Predictions**: Returns average CpG count across entire sequence
- **Track Predictions**: Returns per-base CpG density using sliding window approach
- **Flexible Scaling**: Linear or log2-scaled outputs (default is linear)
- **Multiple Formats**: Supports JSON and MessagePack request/response formats
- **Containerized**:Apptainer container for reproducible deployment

## Quick Start

[To run pre-built CpG Predictor download here](https://zenodo.org/records/18182725)

### Run CpG Prediction
```bash
apptainer run --containall run cpg_predictor.sif <ip_address> <port>
# Example: apptainer run --containall cpg_predictor.sif 0.0.0.0 5000
```

### Build Container
```bash
apptainer build cpg_predictor.sif predictor.def
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
  "readout": "point" | "track",
  "prediction_tasks": [
    {
      "name": "task_name",
      "type": "experiment_type",
      "cell_type": "cell_type",
      "species": "species",
      "scale": "linear" | "log"
    }
  ],
  "sequences": {
    "seq_id": "ACGTACGT..."
  },
  "prediction_ranges": {
    "seq_id": [start, end]
  }
}
```

**Response Format:**
```json
{
  "predictor_name": "CpG Predictor",
  "bin_size": 1,
  "prediction_tasks": [
    {
      "name": "task_name",
      "scale_prediction_actual": "linear",
      "type_requested": "accessibility",
        "type_actual": [
            "NA"
        ],
        "cell_type_requested": "cell_type",
        "cell_type_actual": "NA",
        "species_requested": "mus_musculus",
        "species_actual": "NA",
        "scale_prediction_requested": "linear",
        "scale_prediction_actual": "linear",
        "predictions": {
            "seq_id": [0.025, 0.032, ...]
      }
    }
  ]
}
```

## Algorithm Details

- **Point Mode**: Calculates mean CpG frequency across entire sequence
- **Track Mode**: Uses 50bp sliding window to compute per-base CpG density
- **CpG Definition**: Dinucleotide sequence "CG"
- **Density Calculation**: (CpG count / window length) × 100
- **Prediction ranges**: Sequences are trimmed to the start and end of the ranges before calculating CpG content


## Project Structure

```
CpG_Predictor/
├── predictor_RestAPI.py       # Main Flask API server
├── cpg_utils.py                # CpG calculation functions
├── predictor_help_message.json # Metadata
├── error_checking_functions.py # Validation utilities
├── schema_validation.py        # Request schema validation
└── predictor_content_handler.py # Format encoding/decoding
```