echo "Waiting for 5 seconds..."
#sleep 3


echo "Connect to localhost" >> log.log 
ssh -N -J ybhatti@int4-pub.snellius.surf.nl ybhatti@fcn73.local.snellius.surf.nl -L 51525:localhost:51525 &
echo "connected to local host ${USER}@${LOGIN_HOST} ${USER}@${BATCH_HOST} -L ${PORT}:localhost:${PORT}"

