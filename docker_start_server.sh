#!/bin/bash

cd server

echo "You will enter a screen. Use CTRL+A then CTRL+D to detach screen."

if [ ! $(screen -list | grep -q "server") ] 
then
  screen -S server -dm bash -c "npm start"
fi

screen -R server
