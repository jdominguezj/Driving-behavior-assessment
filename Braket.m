%frame deteciron
clear, clc
tic
f = VideoReader('2019-02-14 15-55-12.flv');
vr = [260.5 648.5 3 34];
va = [274.5 648.5 3 34];
%%
start=14*30;
finish=164*30;
t=start/30:1/30:finish/30;
%72 sec
brake = zeros(1,length(start:finish));
n = 0;

% figure 
for iFrame = start:finish
    n = n+1;
    frame = rgb2gray(read(f, iFrame));
    framec = imcrop(frame, vr);
    framea = imcrop(frame, va);
    frameb = framec>80;
    frameac= framea>80;
    brake(n) = sum(frameb(:))/140;
    acce(n) =sum(frameac(:))/140;
end

%%

figure,hold on,plot(brake),plot(acce)
legend('Braking','Acceleration')

toc