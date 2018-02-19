cp infiles/Reynier/3H_fast_1day_rad_bound.data infiles/current.data 
./simc
python python/make_root_one_default.py
cd output/
mv current.root fast_H_bound.root
cd ..
cp infiles/Reynier/3H_fast_1day_rad_cont.data infiles/current.data 
./simc
python python/make_root_one_default.py
cd output/
mv current.root fast_H_continuum.root
cd ..
cp infiles/Reynier/3H_slow_8day_rad_bound.data infiles/current.data
./simc
python python/make_root_one_default.py
cd output/
mv current.root slow_H_bound.root 
cd ..
cp infiles/Reynier/3H_slow_8day_rad_cont.data infiles/current.data
./simc
python python/make_root_one_default.py
cd output/
mv current.root slow_H_continuum.root
