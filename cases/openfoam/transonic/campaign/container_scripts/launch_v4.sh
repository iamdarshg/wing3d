while ps -ef | grep -v grep | grep -q "solver4[.]sh 1 "; do sleep 45; done
bash /tmp/solver_v4.sh
