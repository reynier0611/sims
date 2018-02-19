cp infiles/Reynier/3He_fast_1day_rad_benhar_bound.data infiles/current.data
./simc
python python/make_root_one_default.py
cd output/
mv current.root fast_He_benhar_bound.root 
cd ..
cp infiles/Reynier/3He_fast_1day_rad_benhar_cont.data infiles/current.data
./simc
python python/make_root_one_default.py
cd output/
mv current.root fast_He_benhar_continuum.root
cd ..
