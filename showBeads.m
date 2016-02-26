function showBeads(cfg)

%cfg=getStrainEnergyAssayConfig();
bp=dlmread(cfg.startposfilename,'\t');
bp=bp(:, 3:5);

t_idx = cfg.tstart;

for z_idx = 0:1000
    if(~exist(cfg.FLfileName(z_idx,t_idx),'file')) 
        break;
    end
        
    img = imread(cfg.FLfileName(z_idx, t_idx));
    img = imadjust(img);

    % current z
    sel = bp(bp(:, 3)==z_idx, :);
    s = size(sel);
    n_beads = s(1);
    
    img = insertShape(img, 'circle', [sel(:, 1:2) 6*ones(n_beads, 1)], 'LineWidth', 1, 'Color', 'yellow');

    % neighbor z
    sel = bp(bp(:, 3)==z_idx-1|bp(:, 3)==z_idx+1, :);
    s = size(sel);
    n_beads = s(1);
    
    img = insertShape(img, 'circle', [sel(:, 1:2) 6*ones(n_beads, 1)], 'LineWidth', 1, 'Color', 'green');
    

    imwrite(img, sprintf('beads_%d_%d.png', t_idx, z_idx));
end

end

