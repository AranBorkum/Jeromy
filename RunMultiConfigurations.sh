for i in 1
do	 
    for j in 1 10 20 30 40 50 60 70 80 90 100 200 300 400 500
    do
	./Jeromy/scripts/generate_likelihood_profiles_cno_discovery.py -s optimistic_uncertainties_hz -m high -e $i -br $j
    done
done
