echo "# " > qwik.cmd 
echo "close all" >> qwik.cmd 
echo "close session" >> qwik.cmd
echo "# " >> qwik.cmd 

#ls *t-10*c-50_*b-205*o-10/clust-0-*trial-9*n-*opt*.pdb | grep -v 'xp' | grep -v 'xm' | sed 's/out/open out/g' >> qwik.cmd
#ls *t-10*c-2*b-1_l-5*/clust-0-*trial-9*n-51*opt*.pdb | grep -v 'xp' | grep -v 'xm' | sed 's/out/open out/g' >> qwik.cmd

ls *t-30*/clust-0-*trial-29*n-*opt*.pdb | grep -v 'xp' | grep -v 'xm' | sed 's/out/open out/g' >> qwik.cmd
#ls *t-10*c-5_*b-5*/clust-1-*trial-9*.pdb | grep -v 'opt'| grep -v 'xp' | grep -v 'xm' | sed 's/out/open out/g' >> qwik.cmd


echo "#">> qwik.cmd
echo "select #0" >> qwik.cmd 
echo "#rlabel sel">> qwik.cmd
echo "runscript showoutline.py #0">> qwik.cmd
echo "#runscript showoutline.py #1">> qwik.cmd
echo "open settings.cmd">> qwik.cmd
echo "turn x -90">>qwik.cmd
echo "turn y -90">>qwik.cmd
echo "focus ;~sel" >> qwik.cmd
echo "#open assigntriangles.center.cmd">>qwik.cmd
echo "#open assigntriangles.cmd">>qwik.cmd
echo "#">> qwik.cmd
echo "#">> qwik.cmd
echo "#">> qwik.cmd




