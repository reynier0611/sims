cp infiles/Igor_d_kin1.data infiles/current.data
./simc
python python/make_root_one_default.py
mv output/current.root output/Igor_d_kin1.root
cp infiles/Igor_3He_bound_kin1.data infiles/current.data
./simc
python python/make_root_one_default.py
mv output/current.root output/Igor_3He_bound_kin1.root
cp infiles/Igor_3He_cont_kin1.data infiles/current.data
./simc
python python/make_root_one_default.py
mv output/current.root output/Igor_3He_cont_kin1.root
