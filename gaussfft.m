clear all; close all; % clear all variables and figures
L=20;	% define the computational domain [-L/2,L/2] 
n=128; % define the number of Fourier modes 2^n

x2=linspace(-L/2,L/2,n+1); % define the domain discretization 

x=x2(1:n); % consider only the first n points: periodicity

u=exp(-x.*x);	% function to take a derivative of 

ut=fft(u);	% FFT the function 

utshift=fftshift(ut); % shift FFT

figure(1), plot(x,u)	% plot initial gaussian 
figure(2), plot(abs(ut))	% plot unshifted transform 
figure(3), plot(abs(utshift)) 