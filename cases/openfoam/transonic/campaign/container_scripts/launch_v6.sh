while ps -ef | grep -v grep | grep -q "solver4[.]sh 0 "; do sleep 45; done
bash /tmp/solver_v6.sh
