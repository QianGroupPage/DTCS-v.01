from surface_crns import SurfaceCRNQueueSimulator

def main():
	# Runs the "Rule 110" example.
	manifest_filename = "rule_110.txt"
	SurfaceCRNQueueSimulator.simulate_surface_crn(manifest_filename)

if __name__ == "__main__":
	main()