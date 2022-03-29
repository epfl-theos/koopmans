from koopmans import io

# Load the workflow object
wf = io.read('ozone.kwf')

# Access the results from the very last calculation
results = wf.calculations[-1].results

# Calculate the IP and EA
ip = -results['homo_energy']
ea = -results['lumo_energy']

# Print
print(f' IP = {ip:.2f} eV')
print(f' EA = {ea:.2f} eV')
