cp infiles/Reynier/3H_slow_1day_rad_kaptari_cont.data infiles/current.data
./simc
python python/make_root_one_default.py
cd output/
mv current.root slow_H_kaptari_cont.root 
cd ..
