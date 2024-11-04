## Files and Structure

- **`gwas_lib.py`**  
  Contains all essential functions.
  
- **`client.py`**  
  This file represents the client-side code and is run multiple times to simulate a multi-client environment. Each client instance performs computations independently, interacting with the server. The script is executed multiple times with multiple calls of `run_client.sh`.

- **`server.py`**  
  This file handles the server-side logic, receiving and aggregating information from each client. It manages secure data processing and combines the results from all clients. The server is activated using `run_server.sh`.

## Bash Scripts

- **`run_client.sh`**  
  A script to initialize each instance of `client.py`, simulating different client environments. Each client is treated as an independent site in the multi-site GWAS study.

- **`run_server.sh`**  
  A script that starts the `server.py`.

- **`server_client.sh`**  
  This overarching bash file is used to initiate both `run_server.sh` and multiple instances of `run_client.sh`.

---
