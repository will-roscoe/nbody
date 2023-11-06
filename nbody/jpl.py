import requests
import re
inp = input('object')
req="https://ssd.jpl.nasa.gov/api/horizons.api?format=json&COMMAND=%27499%27&OBJ_DATA=%27YES%27&MAKE_EPHEM=%27YES%27&EPHEM_TYPE=%27OBSERVER%27&CENTER=%27500@399%27&START_TIME=%272006-01-01%27&STOP_TIME=%272006-01-20%27&STEP_SIZE=%271%20d%27&QUANTITIES=%271,9,20,23,24,29%27"
r = requests.get(req)
#print(f"json={r.json()}")
jsonResponse=r.json()

#for key, value in jsonResponse.items():
#        print(key, ":", value)
        
result = jsonResponse['result']
print(result)
m = re.search(r"\s+?Mass\s+?([\S]+).*?=\s+?([\S]+?)\s+?", result)
r = re.search(r'\s+?radius\s+?([\S]+).*?=\s+?([\S]+?)\s+?', result)
if m:  
    print("m=", m)
    mbase, mexp_s = m.groups()[1], m.groups()[0]
    print("r=", r)
    rad =r.groups()[1]
else:
    print("Mass not found")
print()
mass = float(mbase) * 10**(float(mexp_s.split('x10^')[1]))
radius = float(rad.split('+-')[0])*1000

print(mass, radius)