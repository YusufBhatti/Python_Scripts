echo "Waiting for 5 seconds..."
sleep 5


echo "Connect to localhost" 
ssh -N -J ybhatti@int5-pub.snellius.surf.nl ybhatti@tcn718.local.snellius.surf.nl -L 51525:localhost:51525 ssh -N -J ${USER}@${LOGIN_HOST} ${USER}@${BATCH_HOST} -L ${PORT}:localhost:${PORT} &
echo "connected to local host ${USER}@${LOGIN_HOST} ${USER}@${BATCH_HOST} -L ${PORT}:localhost:${PORT}"

