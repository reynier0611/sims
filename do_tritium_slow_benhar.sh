cp infiles/Reynier/3H_slow_1day_rad_benhar_cont.data infiles/current.data
./simc
python python/make_root_one_default.py
cd output/
mv current.root slow_H_benhar_cont.root 
cd ..
