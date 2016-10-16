graphs : output.csv EnergyVals.csv 
	python EMagCheckerboard.py htmlplot=1

output.csv :
	python EMagCheckerboard.py >> output.csv

EnergyVals.csv :
	perl Postprocess.pl

clean : 
	rm -rf _output
	rm -rf _plots
	rm -f output.csv
	rm -f EnergyVals.csv
