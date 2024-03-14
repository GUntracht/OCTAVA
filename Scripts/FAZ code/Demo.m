% Call the Chan-Vese method in the following way:
% phi = chanvese(I,n,miu,lambda1,lambda2,maxsize)
% Where:
% I is the Image 
% n is the number of iterations
% miu, lambda1 and lambda2 are parameters for Chan-Vese
% maxsize is the max number of pixels on the long edge on the output image

clc;
clear;
close all;
%image should be RGB 
I = imread('BrGe_OD.png');
I = rgb2gray(I);
I = 256 - I;
[phi,s] = chanvese(I,5500,20,5,1,500);
%Input image,iteration times,parameter1,parameter2,parameter3,max size

