%% Channel Repair
% Find neighbours
cfg_neighbour.method = 'template';
cfg_neighbour.layout = 'biosemi64.lay';
cfg_neighbour.channel = repair_electrodes;

cfg_repair.neighbours = ft_prepare_neighbours(cfg_neighbour);

% Establish 3D Electrode Layout
cfg_repair.elec = ft_read_sens('standard_1020.elc'); % this ships with fieldtrip

% Remove electrodes not in 64-elec cap
rmelecs = ~ismember(cfg_repair.elec.label, prep_data.label);
cfg_repair.elec.chanpos(rmelecs, :) = [];
cfg_repair.elec.elecpos(rmelecs, :) = [];
cfg_repair.elec.label(rmelecs) = [];

% "Repair" channels
cfg_repair.method = 'nearest';
cfg_repair.badchannel = repair_electrodes;

interp_data = ft_channelrepair(cfg_repair, prep_data);