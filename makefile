graphs :
	python EMagCheckerboard.py htmlplot=1 >> output.csv
	perl Postprocess.pl
	perl Postprocess2.pl
	cp EnergyVals.csv ./EnergyGraphs/EVal1.csv

clean : 
	rm -rf _output
	rm -rf _plots
	rm -f output.csv
	rm -f LimitCurve.csv
	rm -f EnergyVals.csv
