cp infiles/Reynier/d_fast_1day_rad.data infiles/current.data
./simc
python python/make_root_one_default.py
cd output/
mv current.root d_fast_1day_rad.root 
cd ..
cp infiles/Reynier/d_slow_1day_rad.data infiles/current.data
./simc
python python/make_root_one_default.py
cd output/
mv current.root d_slow_1day_rad.root   
cd ..
