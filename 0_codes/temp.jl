y=zeros(250001,131);
fl="bipolar_pre"
for id in 0:131
	t, y[:,id+1] = read_channel(id, 100, 200.0, fl)    # time series starts at 4e-4
end
	
t, aux = read_channel(0, 100, 200.0, fl)
