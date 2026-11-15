
import time
from genomics.ncbi_annotation import fetch_ncbi_record

start = time.time()

print("Connecting to NCBI...", flush=True)

try:
    print("Retrieving CP002686.1...", flush=True)

    record = fetch_ncbi_record("CP002686.1")

    print("Retrieval successful!", flush=True)
    print("Accession:", record.id)
    print("Sequence length:", len(record.seq))
    print("Features:", len(record.features))
    print("Time taken:", round(time.time() - start, 2), "seconds")

except Exception as error:
    print("TEST FAILED:", error)
