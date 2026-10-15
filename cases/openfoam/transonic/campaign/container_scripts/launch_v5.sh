while ps -ef | grep -v grep | grep -q "solver5[.]sh"; do sleep 45; done
bash /tmp/solver_v5.sh
