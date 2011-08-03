function status=martin_corr_gui4(input,tfmselect,status)

% %version MartinSchorb 100119
% %
% %usage is martin_corr_gui(struct from martin_ls_blind2); status is needed
% as status variable
% %
%
% graphical interface to select transformation to correlate images


if tfmselect>0


%generate image list
a=input.blind(1,tfmselect).image;
% % num=double(imread('/struct/briggs/schorb/code/1.tif'));
% num=double(imread(['/struct/briggs/schorb/code/',int2str(tfmselect),'.tif']));
% num=uint8(num);
% num(:,:,1)=num;
% % num(:,:,2)=num(:,:,1);
% num(:,:,3)=num(:,:,1);
% num(2048,2048,:)=00;
% err=figure('Visible','off','units','Pixels','Position',[0 0 400 200],'units','Pixels','position',[0 0 800 200]);
% axis off;
% as=text(0,0,['closest bead acc.: ',num2str(input.blind(1).closerr),' pixels'],'FontSize',64, 'Color','red');
% errim=getframe(err);
% errim2=errim.cdata-188;
% errim2=errim2/max(max(max(errim2)))*255;
% errim2(2048,2048,:)=0;
% errim2(1:220,1:450,:)=errim2(51:270,101:550,:);
% close(err);
% 
% num(2048,2048,:)=0;
% num=(num+errim2)/2;
% %     num(:,:,1)=imadjust(num(:,:,1));
% 
% a(1:160,1:430,:)=a(1:160,1:430,:)/2+num(1:160,1:430,:);
% k=uint8(188*ones(2048,500,3));
% popup=' 1';
% for ind=2:size(input.blind,2)
%     
%     
%     
%     
%     
%     num=double(imread(['/struct/briggs/schorb/code/',int2str(ind),'.tif']));
%     num=uint8(num);
%     num(:,:,1)=num;
% %     num(:,:,2)=num(:,:,1);
%     num(:,:,3)=num(:,:,1);
%    if isa(input.blind(1,ind).image,'uint8') %find bad transformations
%     %show estimated error
%     err=figure('Visible','off','units','Pixels','Position',[0 0 400 200],'units','Pixels','position',[0 0 800 200]);
%     axis off;
%     as=text(0,0,['Predicted acc.: ',num2str(input.blind(ind).closerr),' pixels'],'FontSize',64, 'Color','red');
%     errim=getframe(err);
%     errim2=errim.cdata-188;
%     errim2=errim2/max(max(max(errim2)))*255;
%     errim2(2048,2048,:)=0;
%     errim2(1:220,1:450,:)=errim2(51:270,101:550,:);
%     close(err);
%         
%     num(2048,2048,:)=0;
%     num=(num+errim2)/2;
% %     num(:,:,1)=imadjust(num(:,:,1));
%     b=input.blind(1,ind).image;
%     b(1:160,1:430,:)=b(1:160,1:430,:)/2+num(1:160,1:430,:);
%     a=[a,k,b];
%     
%     if ind<10
%         popup=[popup;' ',int2str(ind)];
%     else
%         popup=[popup;int2str(ind)];
%     end
%    else
%     num(2048,600,:)=0;
%     a=[a,k,num];
%     popup=[popup;'--'];
%    end
%    
% end

seltext=int2str(tfmselect);

else
% num=double(imread(['/struct/briggs/schorb/code/all.tif']));
%     num=uint8(num);
%     num(:,:,1)=num;
% %     num(:,:,2)=num(:,:,1);
%     num(:,:,3)=num(:,:,1);
%     %show estimated error
%     err=figure('Visible','off','units','Pixels','Position',[0 0 400 200],'units','Pixels','position',[0 0 800 200]);
%     axis off;
%     as=text(0,0,['closest bead acc.: ',num2str(input.all.closerr),' pixels'],'FontSize',64, 'Color','red');
%     errim=getframe(err);
%     errim2=errim.cdata-188;
%     errim2=errim2/max(max(max(errim2)))*255;
%     errim2(2048,2048,:)=0;
%     errim2(1:220,1:450,:)=errim2(51:270,101:550,:);
%     close(err);
    
%     num(:,:,1)=imadjust(num(:,:,1));
    if isa(input.all.image,'uint8')
%         num(2048,2048,:)=0;
%         num=(num+errim2)/2;
        b=input.all.image;
%         b(1:160,1:430,:)=b(1:160,1:430,:)/2+num(1:160,1:430,:);
        a=b;
        seltext='all';
%         popup=[{,popup,'all'}];
    else
       num(2048,600,:)=0;
       a=[];
       seltext='no image';
%        popup=[{,popup,'---'}];
    end
    
end


f = figure('Visible','off','Position','NumberTitle','off',[0,120,900,900]);
s=imshow(a);
hSP = imscrollpanel(f,s);
api = iptgetapi(hSP);
api.setMagnification(0.4);

%  htext = uicontrol('Style','text','String','Select Transformation to use',...
%            'Position',[180,90,450,15]);
%  hpopup = uicontrol('Style','popupmenu',...
%            'String',popup,...
%            'Position',[300,50,100,25],'Callback',{@popup_menu_Callback});

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

