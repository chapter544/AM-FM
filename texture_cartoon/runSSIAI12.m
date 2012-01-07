clear all
close all
clc

inImageDir = '/home/chuong/Images';
outImageDir = '/home/chuong/Research/AMFM/experiments/runChristmas2011/2012_01_05';
M = 256;
N = 256;

% snake1R
cartoonTextureDriver('snake1R', 'float', 256, 256,  inImageDir, outImageDir);

% burlap 
cartoonTextureDriver('burlapR', 'float', 256, 256,  inImageDir, outImageDir);

% house 256x256
cartoonTextureDriver('house', 'png', 0, 0,  inImageDir, outImageDir);

% peppers 256x256
cartoonTextureDriver('peppers256', 'png', 0, 0,  inImageDir, outImageDir);

% lena 512x512 png
cartoonTextureDriver('lena', 'png', 0, 0,  inImageDir, outImageDir);

% boat 512x512 png
cartoonTextureDriver('boat', 'png', 0, 0,  inImageDir, outImageDir);

% fingerprint 512x512
cartoonTextureDriver('fingerprint', 'png', 0, 0,  inImageDir, outImageDir);

% barbara 512x512
cartoonTextureDriver('barbara', 'png', 0, 0,  inImageDir, outImageDir);

% flinttones 512x512
cartoonTextureDriver('flinstones', 'png', 0, 0,  inImageDir, outImageDir);

% kodim01 768x512
cartoonTextureDriver('kodim01', 'png', 0, 0,  inImageDir, outImageDir);

% kodim05
cartoonTextureDriver('kodim05', 'png', 0, 0,  inImageDir, outImageDir);

% kodim08
cartoonTextureDriver('kodim08', 'png', 0, 0,  inImageDir, outImageDir);

% kodim17
cartoonTextureDriver('kodim17', 'png', 0, 0,  inImageDir, outImageDir);

% kodim22
cartoonTextureDriver('kodim22', 'png', 0, 0,  inImageDir, outImageDir);

% kodim23
cartoonTextureDriver('kodim23', 'png', 0, 0,  inImageDir, outImageDir);
