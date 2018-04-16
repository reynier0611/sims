cp infiles/Reynier/3He_slow_1day_rad_benhar_bound.data infiles/current.data
./simc
python python/make_root_one_default.py
cd output/
mv current.root slow_He_benhar_bound.root 
cd ..
cp infiles/Reynier/3He_slow_1day_rad_benhar_cont.data infiles/current.data
./simc
python python/make_root_one_default.py
cd output/
mv current.root slow_He_benhar_continuum.root
