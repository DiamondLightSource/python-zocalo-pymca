#!/bin/bash
tmux new-session \; \
send-keys 'module load activemq/zocdev' C-m \; \
split-window -h \; \
send-keys 'module load dials; dlstbx.status_monitor --test' C-m \; \
select-pane -t 0 \; \
split-window -v \; \
send-keys 'module load dials; dlstbx.service --test -s DLSDispatcher -v' C-m \; \
select-pane -t 2 \; \
split-window -v \; \
send-keys 'module load python/zocalo ; zocalo.service --test -s DLSPyMcaFitter -v' C-m \;
