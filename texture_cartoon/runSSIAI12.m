clear all
close all
clc

inImageDir = '/home/chuong/Images';
outImageDir = '/home/chuong/Research/AMFM/runChristmas2011/2012_12_14';
M = 256;
N = 256;

% snake1R
alpha = 0.001;
beta = 0.003;
cartoonTextureDriver('snake1R', 'float', 256, 256, alpha, beta, inImageDir, outImageDir);

% lena 256x256
alpha = 0.5;
beta = 0.003;
cartoonTextureDriver('lenaR', 'float', 256, 256, alpha, beta, inImageDir, outImageDir);

% lena 512x512 png
alpha = 0.001;
beta = 0.003;
cartoonTextureDriver('lena', 'png', 256, 256, alpha, beta, inImageDir, outImageDir);

% boat 512x512 png
alpha = 0.001;
beta = 0.5;
cartoonTextureDriver('boat', 'png', 256, 256, alpha, beta, inImageDir, outImageDir);

% fingerprint
alpha = -0.5;
beta = 0.5;
cartoonTextureDriver('fingerprint', 'png', 256, 256, alpha, beta, inImageDir, outImageDir);

% house 
alpha = 0.5;
beta = 0.5;
cartoonTextureDriver('house', 'png', 256, 256, alpha, beta, inImageDir, outImageDir);

% peppers 
alpha = 0.001;
beta = 0.003;
cartoonTextureDriver('peppers256', 'png', 256, 256, alpha, beta, inImageDir, outImageDir);

% barbara
alpha = 0.001;
beta = 0.003;
cartoonTextureDriver('barbara', 'png', 256, 256, alpha, beta, inImageDir, outImageDir);

% flinttones
alpha = 0.5;
beta = 0.003;
cartoonTextureDriver('flinstones', 'png', 256, 256, alpha, beta, inImageDir, outImageDir);

% kodim01
alpha = 0.001;
beta = 0.003;
cartoonTextureDriver('kodim01', 'png', 256, 256, alpha, beta, inImageDir, outImageDir);

% kodim05
alpha = 0.001;
beta = 0.003;
cartoonTextureDriver('kodim05', 'png', 256, 256, alpha, beta, inImageDir, outImageDir);

% kodim08
alpha = 0.001;
beta = 0.003;
cartoonTextureDriver('kodim08', 'png', 256, 256, alpha, beta, inImageDir, outImageDir);

% kodim17
alpha = 0.001;
beta = 0.003;
cartoonTextureDriver('kodim17', 'png', 256, 256, alpha, beta, inImageDir, outImageDir);

% kodim22
alpha = 0.001;
beta = 0.003;
cartoonTextureDriver('kodim22', 'png', 256, 256, alpha, beta, inImageDir, outImageDir);

% kodim23
alpha = 0.001;
beta = 0.003;
cartoonTextureDriver('kodim23', 'png', 256, 256, alpha, beta, inImageDir, outImageDir);
