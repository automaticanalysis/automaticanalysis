function stop = pewarp_fminunc_regularised_plot(x, optimValues, state)

global bX bY dbY bZ target8 source mask krn xscale regstrength

[f, H, warped, D] = mex_pewarpcost_regularised(single(x), bX, bY, dbY, bZ, target8, source, mask, krn, xscale, regstrength);

h = findobj('Name', 'pewarpOut');
if ~ishandle(h)
    h = figure('Name', 'pewarpOut', 'Position', [100 500 1200 300]);
end
set(0, 'currentFigure', h);

disp('Drawing...');

colormap('gray');

subplot(2, 3, 1);
image(squeeze(target8(floor(end / 2), :, :))' / 4);
axis xy;

subplot(2, 3, 2);
image(squeeze(warped(floor(end / 2), :, :))' / 4);
axis xy;

subplot(2, 3, 3);
image(squeeze(source(floor(end / 2), :, :))' / 4);
axis xy;

subplot(2, 3, 4);
image(32 + (double(squeeze(warped(floor(end / 2), :, :))') - squeeze(source(floor(end / 2), :, :))') / 4);
axis xy;

subplot(2, 3, 5);
image(squeeze(D(floor(end / 2), :, :))', 'CDataMapping', 'scaled');
axis xy;
set(gca, 'CLim', [-5 5]);
% colorbar;

subplot(2, 3, 6);
image(H, 'CDataMapping', 'scaled');
axis xy;
set(gca, 'CLim', [0 50]);

refresh;
pause(.1);

stop = false;