function aaworker_pollsetup
global aaworker;

aaworker.polltimer=timer('TimerFcn',@aaworker_poll,'Period',2,'ExecutionMode','FixedDelay'); 
start(aaworker.polltimer);

