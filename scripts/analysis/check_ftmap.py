# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
import requests, time

print("Checking FTMap availability...")
try:
    r = requests.get("https://ftmap.bu.edu/", timeout=10)
    print(f"FTMap status: {r.status_code}")
    ct = r.headers.get("content-type", "")
    print(f"Content-Type: {ct}")
    if "form" in r.text.lower():
        print("Has form (web-based submission)")
    if "api" in r.text.lower():
        print("Has API reference")
except Exception as e:
    print(f"FTMap error: {e}")

print("\nChecking proteins.plus for DoGSiteScorer...")
try:
    # Try their documented REST format
    r = requests.post(
        "https://proteins.plus/api/v2/dogsite",
        data={"pdb_code": "6MEO"},
        headers={"Accept": "application/json"},
        timeout=15
    )
    print(f"DoGSiteScorer POST status: {r.status_code}")
    ct = r.headers.get("content-type", "")
    print(f"Content-Type: {ct}")
    if "json" in ct:
        import json
        print(json.dumps(r.json(), indent=2)[:500])
    else:
        # Try as a regular form submit
        r2 = requests.post(
            "https://proteins.plus/api/v2/molecule_handler",
            data={"pdb_code": "6MEO"},
            headers={"Accept": "application/json"},
            timeout=15
        )
        print(f"\nmolecule_handler POST: {r2.status_code}")
        ct2 = r2.headers.get("content-type", "")
        print(f"Content-Type: {ct2}")
        if "json" in ct2:
            print(json.dumps(r2.json(), indent=2)[:500])
except Exception as e:
    print(f"DoGSiteScorer error: {e}")
