#!/bin/bash

(>&2 echo "");
for((i=1; i < 256; i++)); do
	(>&2 echo -ne "\r$((($i * 100)/255))");
	for((k=1; (k+i) <= 256; k++)); do
		./main $i $k 7 $RANDOM; status=$?;
		if [ $status -ne 0 ]; then
			(>&2 echo "ERROR");
			exit;
		fi;
	done;
done
(>&2 echo "");
(>&2 echo "SUCCESS");
