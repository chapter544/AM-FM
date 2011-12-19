clear all
close all
clc

inImageDir = '/home/chuong/Images';
outImageDir = '/home/chuong/Research/AMFM/experiments/runChristmas2011/2012_12_19';
M = 256;
N = 256;

% snake1R
cartoonTextureDriver('snake1R', 'float', 256, 256,  inImageDir, outImageDir);

% lena 256x256
cartoonTextureDriver('lenaR', 'float', 256, 256,  inImageDir, outImageDir);

% lena 512x512 png
cartoonTextureDriver('lena', 'png', 256, 256,  inImageDir, outImageDir);

% boat 512x512 png
cartoonTextureDriver('boat', 'png', 256, 256,  inImageDir, outImageDir);

% fingerprint
cartoonTextureDriver('fingerprint', 'png', 256, 256,  inImageDir, outImageDir);

% house 
cartoonTextureDriver('house', 'png', 256, 256,  inImageDir, outImageDir);

% peppers 
cartoonTextureDriver('peppers256', 'png', 256, 256,  inImageDir, outImageDir);

% barbara
cartoonTextureDriver('barbara', 'png', 256, 256,  inImageDir, outImageDir);

% flinttones
cartoonTextureDriver('flinstones', 'png', 256, 256,  inImageDir, outImageDir);

% kodim01
cartoonTextureDriver('kodim01', 'png', 256, 256,  inImageDir, outImageDir);

% kodim05
cartoonTextureDriver('kodim05', 'png', 256, 256,  inImageDir, outImageDir);

% kodim08
cartoonTextureDriver('kodim08', 'png', 256, 256,  inImageDir, outImageDir);

% kodim17
cartoonTextureDriver('kodim17', 'png', 256, 256,  inImageDir, outImageDir);

% kodim22
cartoonTextureDriver('kodim22', 'png', 256, 256,  inImageDir, outImageDir);

% kodim23
cartoonTextureDriver('kodim23', 'png', 256, 256,  inImageDir, outImageDir);
