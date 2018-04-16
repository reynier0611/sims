cp infiles/Reynier/3He_slow_1day_rad_kaptari_bound.data infiles/current.data
./simc
python python/make_root_one_default.py
cd output/
mv current.root slow_He_kaptari_bound.root 
cd ..
cp infiles/Reynier/3He_slow_1day_rad_kaptari_cont.data infiles/current.data
./simc
python python/make_root_one_default.py
cd output/
mv current.root slow_He_kaptari_continuum.root
