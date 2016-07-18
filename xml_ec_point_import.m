function pos=xml_ec_point_import(file)

st=xml2struct(file);

for i=1:length(st.root.rois.roi)    
    pos(i,1)=str2num(st.root.rois.roi{i}.position.pos_x.Text);
	pos(i,2)=str2num(st.root.rois.roi{i}.position.pos_y.Text);
end