# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
import requests, json, time, sys

PDB_ID = '6MEO'

print('=== DoGSiteScorer via proteins.plus API ===')

# Try the documented API
print('Trying proteins.plus REST endpoints...')

# Check API root
r = requests.get('https://proteins.plus', timeout=15)
print(f'Root status: {r.status_code}, type: {r.headers.get("content-type","")}')

# Try submitting to dogsite endpoint
endpoints = [
    'https://proteins.plus/api/v2/dogsite',
    'https://proteins.plus/api/dogsite',
    'https://proteins.plus/api/v2/molecule_handler',
]

for ep in endpoints:
    try:
        r = requests.post(ep, json={'pdb_code': PDB_ID}, timeout=15)
        print(f'\n{ep}')
        print(f'  Status: {r.status_code}')
        print(f'  Content-Type: {r.headers.get("content-type","")}')
        ct = r.headers.get("content-type","")
        if 'json' in ct:
            print(f'  Body: {json.dumps(r.json(), indent=2)[:500]}')
        else:
            print(f'  Body (first 300): {r.text[:300]}')
    except Exception as e:
        print(f'\n{ep}: {e}')

# Try GET on pdb info
try:
    r = requests.get(f'https://proteins.plus/api/v2/pdb/{PDB_ID}', timeout=15)
    print(f'\nGET pdb/{PDB_ID}: status={r.status_code}')
    ct = r.headers.get("content-type","")
    if 'json' in ct:
        print(json.dumps(r.json(), indent=2)[:500])
    else:
        print(r.text[:300])
except Exception as e:
    print(f'GET pdb error: {e}')
