function status=martin_corr_gui4(input,tfmselect,status)

% % version MartinSchorb 130312
% % Copyright EMBL 2013, All rights reserved
% %
% %usage is martin_corr_gui(struct from martin_ls_blind2); status is needed
% as status variable
% %
%
% graphical interface to select transformation to correlate images


if tfmselect>0


else

    if isa(input.all.image,'uint8')

        b=input.all.image;

        a=b;
        seltext='all';

    else
       num(2048,600,:)=0;
       a=[];
       seltext='no image';

    end
    
end


f = figure('Visible','off','Position',[0,120,900,900],'NumberTitle','off');
s=imshow(a);
hSP = imscrollpanel(f,s);
api = iptgetapi(hSP);
api.setMagnification(0.4);


 hgo = uicontrol('Style','pushbutton','String','GO',...
          'Position',[660,90,70,25],...
          'Callback',{@gobutton_Callback});
 hsel = uicontrol('Style','pushbutton','String','select new fiducials',...
          'Position',[660,50,140,25],...
          'Callback',{@selbutton_Callback});
 selname = uicontrol('Style','text','String',['optimal bead configuration: ',seltext],...
          'Position',[300,90,140,45]);

% Assign the GUI a name to appear in the window title.
set(f,'Name','Check the transformation to use for correlation')


% Move the GUI to the center of the screen.
movegui(f,'center')


% Make the GUI visible.
set(f,'Visible','on')

select=1;

%  Pop-up menu callback. Read the pop-up menu Value property to
%  determine which item is currently displayed and make it the
%  current data. This callback automatically has access to 
%  current_data because this function is nested at a lower level.
   function popup_menu_Callback(hObject,eventdata) 
      % Determine the selected data set.
%       str = get(source, 'String')
    select = get(hObject,'Value');
%       guidata(hObject, handles);
      % Set current data to the selected data set.
      
   end

% go button Callback
function gobutton_Callback(hObject,eventdata) 
status=1;
close,f;
% guidata(hObject, handles);
    return
   end

% select again button callback
function selbutton_Callback(hObject,eventdata) 
close,f;
% guidata(hObject, handles);
    return 
end
waitfor(f);
   
return
end

