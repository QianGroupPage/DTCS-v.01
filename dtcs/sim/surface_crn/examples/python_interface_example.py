from surface_crns import SurfaceCRNQueueSimulator

import dtcs.sim.surface_crn.core


def main():
	# Runs the "Rule 110" example.
	manifest_filename = "rule_110.txt"
	dtcs.sim.surface_crn.core._simulate_surface_crn(manifest_filename)

if __name__ == "__main__":
	main()