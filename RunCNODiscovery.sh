for i in 10 20 30 50 100 200 300 400 500 800 1000
do
    for j in 1
    do 
	./Jeromy/scripts/generate_likelihood_profiles_cno_discovery.py -e $j -b $i -m high -s standard_uncertainties_hz
    done
done    
