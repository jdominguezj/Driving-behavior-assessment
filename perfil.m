clear, clc
%Video files

% v = VideoReader('2019-02-12 11-29-51.flv');
% v = VideoReader('2019-02-13 16-33-19.flv'); %start 10 seconds finish 184
% v = VideoReader('2019-02-14 10-30-48.flv');  %start 88 seconds finish 236
%v = VideoReader('2019-02-14 15-00-22.flv');  %start 33 seconds finihs 94
v = VideoReader('2019-02-14 15-55-12.flv');  %start 14 finish 596
 

%Load variables neccesaries

load rec.mat
load 3positions.mat
load digits.mat

%Start and end of the  video

start=14*30;
finish=596*30;

%perfilv = zeros(1,length(840:6600));
%perfilv = zeros(1,length(300:5520));
perfilv = zeros(1,length(start:finish));
%t = 10:1/30:184;

%Time vector
t=start/30:1/30:finish/30

%Initializing variable
n = 0;

figure
for iFrame=14*30:596*30
    n = n+1;
    
    framec = rgb2gray(read(v, iFrame));
    frame = imbinarize(imcrop(framec,rec));
  
    pos1 = imcrop(frame,a);
    pos2 = imcrop(frame,b);
    pos3 = imcrop(frame,c);
    
    resta = zeros(1,10);
    for i = 1:10
        resta(i) = norm(pos1-digits(:,:,i));
    end
    p1 = find(resta == min(resta))-1;
    
    
    resta = zeros(1,10);
    for i = 1:10
        resta(i) = norm(pos2-digits(:,:,i));
    end
    p2 = find(resta == min(resta))-1;
    
    
    resta = zeros(1,2);
    for i = 1:2
        resta(i) = norm(pos3-digits(:,:,i));
    end
    p3 = find(resta == min(resta))-1;
 
    perfilv(n) = 100*p3+10*p2+p1;
    
    framec = insertText(framec,[180 45],perfilv(n),'FontSize',30);
    imagesc(framec), colormap gray, axis off image
    drawnow
    
%     pause(0.2)
end 

%Time vs Speed
matrix= [t perfilv];

%Plot
figure (2), plot(t-min(t),perfilv)

%Time from zero
rt=t-min(t);
%%
% Filtering 

%Mean moving filter for acceleration and speed
mmv=movmean(perfilv,60);
mm=movmean(diff(perfilv),60);
k=mm;

%Generating acceleration and braking events

f=zeros(1,length(k));
ng=zeros(1,length(k));

for i=1:1:length(k)
%k(i) = (k(i)-k(i-29));
if k(i)>=0
    f(i)=1; %Acc
else 
    ng(i)=1; %Brake
end
end

%Plotting results
figure,hold on, plot(ng),plot(mmv),plot(f),plot(mm)
legend('Bool.Brake','Vel','Bool.Acc','Acc')

%Braking events
fv=[rt(1:end-1);ng]';

%Speed profile
vf=[rt;mmv]';

%Print vectors to send to OpenModelica
form='%10.2f,%10.2f;';
sprintf(form,fv')
sprintf(form,vf')

%%
%Braking bar crop
vr = [260.5 648.5 3 34];

%Creating vector for brake pressure
sumas = zeros(1,length(start:finish));
n = 0;

for iFrame = start:finish
    n = n+1;
    frame = rgb2gray(read(f, iFrame));
    framec = imcrop(frame, vr);
    frameb = framec>80;
    sumas(n) = sum(frameb(:))/140;
end

%Braking pressure along time
figure, plot(t,sumas)
