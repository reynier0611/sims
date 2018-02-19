cp infiles/Reynier/3He_fast_1day_rad_bound.data infiles/current.data 
./simc
python python/make_root_one_default.py
cd output/
mv current.root fast_He_bound.root
cd ..
cp infiles/Reynier/3He_fast_1day_rad_cont.data infiles/current.data 
./simc
python python/make_root_one_default.py
cd output/
mv current.root fast_He_continuum.root
cd ..
cp infiles/Reynier/3He_slow_1day_rad_bound.data infiles/current.data
./simc
python python/make_root_one_default.py
cd output/
mv current.root slow_He_bound.root 
cd ..
cp infiles/Reynier/3He_slow_1day_rad_cont.data infiles/current.data
./simc
python python/make_root_one_default.py
cd output/
mv current.root slow_He_continuum.root
