function WaitAndRetunArray(sec, N)
Array = rand(N)
sleep(sec)
return(Array)
end