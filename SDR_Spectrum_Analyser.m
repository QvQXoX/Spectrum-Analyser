function varargout = SDR_Spectrum_Analyser(varargin)
% SDR_SPECTRUM_ANALYSER MATLAB code for SDR_Spectrum_Analyser.fig
%      SDR_SPECTRUM_ANALYSER, by itself, creates a new SDR_SPECTRUM_ANALYSER or raises the existing
%      singleton*.
%
%      H = SDR_SPECTRUM_ANALYSER returns the handle to a new SDR_SPECTRUM_ANALYSER or the handle to
%      the existing singleton*.
%
%      SDR_SPECTRUM_ANALYSER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SDR_SPECTRUM_ANALYSER.M with the given input arguments.
%
%      SDR_SPECTRUM_ANALYSER('Property','Value',...) creates a new SDR_SPECTRUM_ANALYSER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SDR_Spectrum_Analyser_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SDR_Spectrum_Analyser_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SDR_Spectrum_Analyser

% Last Modified by GUIDE v2.5 08-Mar-2016 01:33:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SDR_Spectrum_Analyser_OpeningFcn, ...
                   'gui_OutputFcn',  @SDR_Spectrum_Analyser_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% EF initialization code - DO NOT EDIT


% --- Executes just before SDR_Spectrum_Analyser is made visible.
function SDR_Spectrum_Analyser_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SDR_Spectrum_Analyser (see VARARGIN)

global off
global dis
global y_max
global y_min
global start_f
global end_f
global centre_f
global f_span
global window
global h
global plot_type 
plot_type = 1;
off = 1;%loop break value, if off == 1, it will stop the loop. 
dis = 1; %Set the default scale mode
y_max = 100;%Set the default max amplitude
y_min = -50;%Set the default min amplitude
start_f = 100;%Set the default start frequency
end_f = 110;%Set the default end frequency
centre_f = 105;%Set the default centre frequency
f_span = 10; %Set the default frequency span
window = 1; %Set the default windowing method

set(handles.SF,'string',start_f);%dispay the default start frequency in edit box
set(handles.EF,'string',end_f);%dispay the default end frequency in edit box
set(handles.CF,'string',centre_f);%dispay the default centre frequency in edit box
set(handles.FS,'string',f_span);%dispay the default frequency span in edit box
set(handles.Y_MAX,'string',y_max);%dispay the default max amplitude in edit box
set(handles.Y_MIN,'string',y_min);%dispay the default min amplitude in edit box
set(handles.FBS,'string','1000');
set(handles.BO,'string','35');

info = sdrinfo();
set(handles.Gain,'string',num2str(info.GainValues));
h = comm.SDRRTLReceiver('CenterFrequency', 100*1e6            ,...
                             'SampleRate',      1e6           ,...
                             'SamplesPerFrame', 256           ,...
                             'EnableTunerAGC',  false         ,...
                             'TunerGain',       10            ,...
                             'OutputDataType',  'double'        ); 

% Choose default command line output for SDR_Spectrum_Analyser
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SDR_Spectrum_Analyser wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SDR_Spectrum_Analyser_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;




% --- Executes during object creation, after setting all properties.
function FS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Start.
function Start_Callback(hObject, eventdata, handles)
% hObject    handle to Start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cla reset % clear the previous figure
global off
global dis
global y_max
global y_min
global data_max
global data_min
global data_avg
global f_span
global start_f
global end_f
global window
global data
global SamplePerFrame
global h
global path
global basename
global sweepperfile
global filetype
global filetype_num
global waterfall_loop
global SweepsPerPlot
global psd
global length_data
global spectrogram_x_axis
global spectrogram_y_axis
global reduction
global plot_type
global mean_time
global noiselevel 
%Disable some parameters---------------------------------------------------

    set(handles.Start,'Enable','Off');
    set(handles.SF,'Enable','Off');             
    set(handles.EF,'Enable','Off');             
    set(handles.CF,'Enable','Off');            
    set(handles.FS,'Enable','Off');              
    set(handles.Auto_Gain,'Enable','Off');      
    set(handles.Gain,'Enable','Off');           
    set(handles.Sample_Frame,'Enable','Off');   
    set(handles.Hardware_Report,'Enable','Off');
    set(handles.FBS,'Enable','Off');
    set(handles.BO,'Enable','Off');
    set(handles.Path,'enable','Off');
    set(handles.BaseName,'enable','Off');
    set(handles.SweepPerFile,'enable','Off');
    set(handles.FileType,'enable','Off');
    set(handles.Browse,'enable','Off');
    set(handles.Save_data,'enable','Off');
    set(handles.Credit,'enable','Off');

%Parameters Setting--------------------------------------------------------

% Frequency Parameters
    contents_SamplePerFrame = cellstr(get(handles.Sample_Frame,'String'));
    SamplePerFrame          = str2double(contents_SamplePerFrame{get(handles.Sample_Frame,'Value')});
    SampleRate              = str2double(get(handles.FBS,'String'))*1000;    % Samples per second
    CenterFrequency         = start_f*1e6 + SampleRate/2;

% File Saving Parameters
    path                    = get(handles.Path,'string');
    basename                = get(handles.BaseName,'string');
    sweepperfile            = str2double(get(handles.SweepPerFile,'string'));
    contents                = cellstr(get(handles.FileType,'String'));
    filetype_num            = get(handles.FileType,'Value');
    filetype                = contents{filetype_num};
    loop                    = 1;
    file_number             = 1;
    
% Ampitude Parameters
    if f_span*1e6 ~= SampleRate 
        % set the lengthe of the data according to the frequency span
        length_data = SamplePerFrame* f_span*1e6/SampleRate;
    else length_data = SamplePerFrame;
    end
    contents_gain           = cellstr(get(handles.Gain,'String'));
    gain                    = str2double(contents_gain{get(handles.Gain,'Value')});
    dis                     = get(handles.display,'Value');
    data_max                = linspace(0,0,length_data); %
    data_min                = linspace(0,0,length_data); %Set default value for three holds.
    data_avg                = linspace(0,0,length_data); %
    scale                   = 1; %Auto scacle factor, When the value of scale changes, 
                                 %adjust the amplitude axis automatically.
    auto_times = 0;
    
% Plot Parameters
    Overlap_factor          = str2double(get(handles.BO,'String'))/100;
    reduction               = ceil(get(handles.Reduc,'Value'));
    SweepsPerPlot           = str2double(get(handles.SweepsPerPlot,'string'));
    waterfall_loop          = 1;
    psd                     = zeros(length_data,SweepsPerPlot);
    spectrogram_x_axis      = zeros(SweepsPerPlot,length_data);
    spectrogram_y_axis      = zeros(SweepsPerPlot,length_data);
    for count=1:SweepsPerPlot 
        % Frequency axis setting for the waterfall plot
        if f_span*1e6 ~= SampleRate
            spectrogram_x_axis(count,:)=linspace(start_f, (start_f*1e6 + ceil((end_f - start_f)*1e6/SampleRate)*SampleRate)/1e6, length_data);
        else
            spectrogram_x_axis(count,:) = linspace(start_f,end_f,length_data);    
        end
    end
    if f_span*1e6 ~= SampleRate
        % block number calculation
        block_num = ceil(f_span*1e6/(SampleRate*(1-Overlap_factor)));
    else block_num = 1; 
    end
    
% Other Parameters
    off                     = 0; % Loop factor, when off = 1, loop breaks.
    noiselevel              = 82;
    loop_number             = 1;% Loop number counter

%--------------------------------------------------------------------------

set(handles.Sample_Rate,'String',num2str(SampleRate/1000));

%DVB-T Paramters setting---------------------------------------------------

h.CenterFrequency           = CenterFrequency;
h.SampleRate                = SampleRate;
h.SamplesPerFrame           = SamplePerFrame;
if get(handles.Auto_Gain,'Value')==1
     h.EnableTunerAGC       = true;
else h.EnableTunerAGC = false;
     h.TunerGain            = gain;
end

%Button Switch-------------------------------------------------------------

%Swith the start button to stop button
drawnow
set(handles.Start,'Visible','Off');
set(handles.Stop,'Visible','On');
set(handles.Start,'Enable','On');

%Main Part-----------------------------------------------------------------

while off==0
      tic
      ph = get(handles.PH,'Value');
      oph = get(handles.OPH,'Value');
      avg = get(handles.AVG,'Value');
  
      %Data capation & processing------------------------------------------     
      
      for i=1:block_num
      
        h.CenterFrequency = start_f*1e6 + (i-1/2)*SampleRate*(1-Overlap_factor);
        [x,~,LOST,LATE] = step(h);
        lost = LOST;
        late = LATE;
        if get(handles.Remove_DC,'Value')==1
           % Remove the DC offset
           x = x - mean(x);
        end
        
        %Windowing Methods-------------------------------------------------
        
        switch window
             case 1  %rectangular
                  x = x.*rectwin(length(x));
             case 2  %Bartlett
                  x = x.*bartlett(length(x));
             case 3  %Bartlett-Hann
                  x = x.*barthannwin(length(x));
             case 4  %Blackman
                  x = x.*blackman(length(x));
             case 5  %Blackman-Harris
                  x = x.*blackmanharris(length(x));
             case 6  %Blackman-Nuttall
                  x = x.*(0.3635819-0.4891775*cos(2*(1:length(x))*pi/length(x))+0.1365995*cos(4*(1:length(x))*pi/length(x))-0.0106411*cos(6*(1:length(x))*pi/length(x)))';
             case 7  %Bohman
                  x = x.*bohmanwin(length(x)); 
             case 8  %Chebyshey
                  x = x.*chebwin(length(x)); 
             case 9  %Cosine
                  x = x.*(sin(pi*(1:length(x))/length(x)))';
             case 10 %Flat top
                  x = x.*flattopwin(length(x));
             case 11 %Gaussian
                  x = x.*gausswin(length(x));
             case 12 %Hamming
                  x = x.*hamming(length(x));
             case 13 %Hann (Hanning)
                  x = x.*hann(length(x));
             case 14 %Kaiser
                  x = x.*kaiser(length(x));
             case 15 %Nuttall
                  x = x.*nuttallwin(length(x));
             case 16 %Parzen
                  x = x.*parzenwin(length(x));
             case 17 %Taylor
                  x = x.*taylorwin(length(x));
             case 18 %Triangular Windowing
                  x = x.*triang(length(x));
             case 19 %Tukey
                  x = x.*tukeywin(length(x));
             case 20 %Welch 
                  x = x.*(1-(2*(1:length(x))/length(x)-1).^2)';  
        end
        
        %linear or logarithmic---------------------------------------------
        
        temp = 20*log10(abs(fftshift(fft(x))))-noiselevel;
        if dis ==1
            temp =10.^(temp/20);
        end
        
        %Smooth function---------------------------------------------------

        if get(handles.Smooth,'Value')==1
            temp = sgolayfilt(temp,3,81);
        end
        
        %Overlap factor----------------------------------------------------
        
        if Overlap_factor~=0
        temp = temp(ceil(SamplePerFrame*Overlap_factor/2):ceil(SamplePerFrame*(1-Overlap_factor/2)));
        end
        if i ==1
            data = temp;
        else data = [data;temp];
        end
        
      end
      
      data = data(1:length_data);
      
      %Max Hold, Average Hold & Min Hold-----------------------------------     
      
      %Max Hold
      if ph==1
         if mean(data_max) == 0
            data_max = data;
                  %Hold the first value until the next loop.
         else data_max = max(data_max, data);
         end
      else data_max = linspace(0,0,length_data);
           %set the sata_ph bake to the default value.
      end
      
      %Min Hold    
      if oph==1
         if mean(data_min) == 0
            data_min = data;
         else data_min = min(data_min, data);
         end
      else data_min = linspace(0,0,length_data);
      end
      
      %Average Hold    
      if avg==1
         if mean(data_avg) == 0
            data_avg = data;
            
            %Hold the first value
            
            avg_n = loop_number;
            %recode the times of loop, when the check box of
            %average hold is clicked.
         else 
               data_avg = data_avg.*((loop_number-avg_n)/(loop_number-avg_n+1)) + data./(loop_number-avg_n+1); 
         end
      else data_avg = linspace(0,0,length_data);
      end
      
      %Power Spectrum Density------------------------------------------------------------------      
      
      x_axis=linspace(start_f, (h.CenterFrequency+SampleRate/2)/1e6, length_data);
      if get(handles.PSD,'value')==1    
         drawnow
         p = plot(handles.axes1,x_axis,data,x_axis,data_avg,x_axis,data_min,x_axis,data_max);
         p(1).Color = [0.3922,0.3922,0.3922];
         p(1).LineWidth = 1;
         if mean(data_avg) == 0
            p(2).LineStyle = 'none';
         %when the average hold is not choosed, set the color as white to
         %make it invisible.
         else p(2).Color = [0.2039,0.6667,0.8627];
         end
         p(2).LineWidth = 1;
         if mean(data_min) == 0
            p(3).LineStyle = 'none';
         else p(3).Color = [0.2980,0.8510,0.3922];
         end
         p(3).LineWidth = 1;
         if mean(data_max) == 0
            p(4).LineStyle = 'none';
         else p(4).Color = [1,0.2314,0.1882];
         end
         p(4).LineWidth = 1;       
         grid on;
         if dis == 1
            ylabel('Amplitude (linear)');
         else ylabel('Amplitude (dB)');
         end
         xlabel('Frequency (MHz)');
         axis([start_f, end_f, y_min, y_max]);
      end 
      %2D & 3D Waterfall-----------------------------------------------------------------------      
      
      for count=1:SweepsPerPlot
          % Time axis setting for waterfall plot 
        if loop_number ==1
            spectrogram_y_axis(count,:) = linspace(count*0.5,count*0.5,length_data);
        else
            spectrogram_y_axis(count,:) = linspace(count*mean_time,count*mean_time,length_data);
        end
      end
      
      %2D Waterfall
      if (get(handles.TwoD,'value')==1)||(get(handles.ThreeD,'value')==1)
          psd(:,(2:SweepsPerPlot)) = psd(:,(1:SweepsPerPlot-1));
          psd(:,1) = data;
          max_psd = max(max(psd));
          min_psd = min(min(psd));
          drawnow
          
%           fig = mesh(handles.axes1,spectrogram_x_axis(:,1:reduction:end),spectrogram_y_axis(:,1:reduction:end),(psd(1:reduction:end,:))','FaceColor','interp',...
%           'EdgeColor','none','EdgeAlpha','0.3',...
%           'FaceLighting','gouraud');
          fig = mesh(spectrogram_x_axis(:,1:reduction:end),spectrogram_y_axis(:,1:reduction:end),(psd(1:reduction:end,:))','FaceColor','interp',...
          'EdgeColor','none','EdgeAlpha','0.3',...
          'FaceLighting','gouraud');
          
          xlabel('Frequency (MHz)');
          ylabel('Time (s)');
          zlabel('Amplitude (Linear)');
          if get(handles.TwoD,'value')==1
             view(0,90);

          else 
              fig.EdgeColor = 'white';
              view(3);
          end
           xlim([start_f end_f]);
           zlim([y_min y_max]);
          
         if loop_number==1
         ylim([0.5 (SweepsPerPlot*0.5)]);
         else ylim([mean_time (SweepsPerPlot*mean_time)]);
         end
            %zlim([y_min y_max]);
          
          waterfall_loop = waterfall_loop + 1;
      end
      %Auto Sacle------------------------------------------------------------------------------      
      %when the display type, windowing mentod, or plot type, any of them
      %changes, it will do auto scale function twice to make sure get the
      %right scale.
      plot_type = get(handles.ThreeD,'value')*4 + get(handles.TwoD,'value')*2 + get(handles.PSD,'value');
      if scale ~= window*10 + dis + plot_type*0.1
          % when scale changes do auto scale.
          scale = 0;
          auto_times = auto_times + 1;
         switch plot_type
             case 1
                 
                if get(handles.PH,'Value')==1
                   y_max = max(data_max);
                else y_max = 1.05*max(max(data))-0.05*min(min(data));
                end
                if get(handles.OPH,'Value')==1
                   y_min = min(data_min);
                else y_min = 1.05*min(min(data))-0.05*max(max(data));
                end
             case 2
                 
                 y_max = 1.05*max_psd-0.05*min_psd;
                 y_min = 1.05*min_psd-0.05*max_psd;
             case 4
                 
                 if loop_number > 0
                     % when the first loop, max value and min value are both 0, 
                     % so i decide to do auto scale at the second loop.
                    if get(handles.display,'Value')==1
                        y_max = 1.05*max_psd-0.05*min_psd;
                        y_min = 0;
                    else
                        y_max = 1.05*max_psd-0.05*min_psd;
                        y_min = -100;
                    end
                 end
         end
         set(handles.Y_MAX,'string',y_max);
         set(handles.Y_MIN,'string',y_min);
      end
      if auto_times >= 2
          % when set the scale twice, stop the auto scale.
          scale = window*10 + dis + plot_type*0.1;
          auto_times = 0;
      end
      
      %Time parameters Dispaly-----------------------------------------------------------------      
      
      loop_time(loop_number) = toc;%recode the time of evey sweep costs
      set(handles.Lost,'string',num2str(lost));% display lost
      set(handles.Late,'string',num2str(late));% display latency
      set(handles.Real_time,'string',num2str(loop_time(loop_number)));
      % display the time that this loop costs
      
      if loop_number > 1
          % calculate the everage loop time
          mean_time = mean(loop_time(2:loop_number));
      else mean_time = mean(loop_time);
      end
      set(handles.AVG_Time,'string',num2str(mean_time));
      %display the average loop time
      set(handles.Sweep_num,'string',num2str(loop_number));
      %display the loop number
      
      %File Saving-----------------------------------------------------------------------------      
      
      if get(handles.Save_data,'Value')==1
          if loop <=sweepperfile
              amplitude(:,loop) = data;
              time(:,loop) = clock;
              frequency = x_axis;
              loop = loop +1;
          else
              loop =1;
              amplitude = amplitude';
              time = time';
              
              % define the file saving path
              amplitudename = [path,'/',basename,sprintf('%03d',file_number),'-Amplitude.',filetype];
              timename = [path,'/',basename,sprintf('%03d',file_number),'-Time.',filetype];
              frequencyname = [path,'/',basename,'-Frequency.',filetype];
              
              switch filetype_num
                  case 1
                      %save as mat format
                      save(amplitudename,'amplitude','-ascii');
                      save(timename,'time','-ascii');
                      save(frequencyname,'frequency','-ascii');
                  case 2
                      %save as txt format
                      save(amplitudename,'amplitude');
                      save(timename,'time');
                      save(frequencyname,'frequency');
              end
              
              file_number = file_number+1;
              clear amplitude;
              clear time;
              amplitude(:,loop) = data;
              time(:,loop) = clock;
              loop = loop +1;
          end
      else
          loop = 1;
          file_number = 1;
      end
      if get(handles.ThreeD,'value')==1
          clear data;
          data = psd;
          
          
      end

      %----------------------------------------------------------------------------------------      
      refresh(gcf)    
      loop_number=loop_number+1;
          
          
end
%----------------------------------------------------------------------------------------      

if get(handles.Save_data,'Value')==1
   amplitude = amplitude';
   time = time';
   frequency = x_axis;
   amplitudename = [path,'/',basename,sprintf('%03d',file_number),'-Amplitude.',filetype];
   timename = [path,'/',basename,sprintf('%03d',file_number),'-Time.',filetype];
   frequencyname = [path,'/',basename,'-Frequency.',filetype];
   switch filetype_num
          case 1
               save(amplitudename,'amplitude','-ascii');
               save(timename,'time','-ascii');
               save(frequencyname,'frequency','-ascii');
          case 2
               save(amplitudename,'amplitude');
               save(timename,'time');
               save(frequencyname,'frequency');
   end
   clear amplitude;
   clear time;
end
%----------------------------------------------------------------------------------------      

 release(h);
 
 
 % --- Executes on button press in Stop.
function Stop_Callback(hObject, eventdata, handles)
% hObject    handle to Stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global off;
off = 1;
set(handles.Start,'Visible','On');
set(handles.Stop,'Visible','Off');
set(handles.SF,'Enable','On');
set(handles.EF,'Enable','On');
set(handles.CF,'Enable','On');
set(handles.FS,'Enable','On');
set(handles.Auto_Gain,'Enable','On');
set(handles.Sample_Frame,'Enable','On');
set(handles.Hardware_Report,'Enable','On');
set(handles.Credit,'Enable','On');
if str2double(get(handles.FS,'string'))*1e3~=str2double(get(handles.FBS,'string'))
   set(handles.FBS,'Enable','On');
   set(handles.BO,'Enable','On');
end
if get(handles.Auto_Gain,'Value')==0
    set(handles.Gain,'Enable','On');
end
set(handles.Save_data,'enable','On');
if get(handles.Save_data,'Value')==1
    set(handles.Path,'enable','On');
    set(handles.BaseName,'enable','On');
    set(handles.SweepPerFile,'enable','On');
    set(handles.FileType,'enable','On');
    set(handles.Browse,'enable','On');
end
if get(handles.PSD,'Value')==1
    set(handles.SweepsPerPlot,'Enable','Off');
    set(handles.PH,'Enable','On');
    set(handles.OPH,'Enable','On');
    set(handles.AVG,'Enable','On');
else set(handles.SweepsPerPlot,'Enable','On')
end





function SP_Callback(hObject, eventdata, handles)
% hObject    handle to SP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SP as text
%        str2double(get(hObject,'String')) returns contents of SP as a double


% --- Executes during object creation, after setting all properties.
function SP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function display_Callback(hObject, eventdata, handles)
global dis
global data_max
global data_min
global data_avg
global waterfall_loop
global psd
global length_data
global SweepsPerPlot
waterfall_loop = 1;
if get(handles.ThreeD,'value')==1
    psd=(-100).*ones(length_data,SweepsPerPlot);
else
    psd=zeros(length_data,SweepsPerPlot);
end
data_max = linspace(0,0,length_data); 
data_min = linspace(0,0,length_data);
data_avg = linspace(0,0,length_data);
dis = get(hObject,'Value');


% --- Executes during object creation, after setting all properties.
function display_CreateFcn(hObject, eventdata, handles)
% hObject    handle to display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PH.
function PH_Callback(hObject, eventdata, handles)
% hObject    handle to PH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of PH


% --- Executes on selection change in popupmenu4.

function AVG_Callback(hObject, eventdata, handles)
% hObject    handle to AVG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of AVG


% --- Executes on button press in OPH.
function OPH_Callback(hObject, eventdata, handles)
% hObject    handle to OPH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of OPH



function Y_MAX_Callback(hObject, eventdata, handles)
% hObject    handle to Y_MAX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global y_max
y_max = str2double(get(hObject,'String'));
% Hints: get(hObject,'String') returns contents of Y_MAX as text
%        str2double(get(hObject,'String')) returns contents of Y_MAX as a double


% --- Executes during object creation, after setting all properties.
function Y_MAX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Y_MAX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Y_MIN_Callback(hObject, eventdata, handles)
% hObject    handle to Y_MIN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global y_min
y_min = str2double(get(hObject,'String'));
% Hints: get(hObject,'String') returns contents of Y_MIN as text
%        str2double(get(hObject,'String')) returns contents of Y_MIN as a double


% --- Executes during object creation, after setting all properties.
function Y_MIN_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Y_MIN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in Auto_Scale.
function Auto_Scale_Callback(hObject, eventdata, handles)
% hObject    handle to Auto_Scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global data_max
global data_min
global data
global y_max
global y_min
global psd
global plot_type
global noiselevel
if get(handles.ThreeD,'value')==1
    clear data;
    data = psd;
end
switch plot_type
      case 1
          
          if get(handles.PH,'Value')==1
             y_max = max(data_max);
          else y_max = 1.05*max(max(data))-0.05*min(min(data));
          end
          if get(handles.OPH,'Value')==1
             y_min = min(data_min);
          else y_min = 1.05*min(min(data))-0.05*max(max(data));
          end
       case 2
           
           y_max = 1.05*max(max(psd))-0.05*min(min(psd));
           y_min = 1.05*min(min(psd))-0.05*max(max(psd));
       case 4
           
           
           if get(handles.display,'Value')==1
               y_max = 1.05*max(max(psd))-0.05*min(min(psd));
               y_min = 0;
           else
               y_max = 1.05*max(max(psd))-0.05*min(min(psd));
               y_min = -noiselevel;
           end
end
set(handles.Y_MAX,'string',y_max);
set(handles.Y_MIN,'string',y_min);
clear data;


% --- Executes on selection change in Windowing.
function Windowing_Callback(hObject, eventdata, handles)
% hObject    handle to Windowing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global window
global data_max
global data_min
global data_avg
global waterfall_loop
global psd
global length_data
global SweepsPerPlot
waterfall_loop = 1;
if get(handles.ThreeD,'value')==1
    psd=(-100).*ones(length_data,SweepsPerPlot);
else
    psd=zeros(length_data,SweepsPerPlot);
end
data_max = linspace(0,0,length_data); 
data_min = linspace(0,0,length_data);
data_avg = linspace(0,0,length_data);
window = get(hObject,'Value');


% Hints: contents = cellstr(get(hObject,'String')) returns Windowing contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Windowing


% --- Executes during object creation, after setting all properties.
function Windowing_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Windowing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function SF_Callback(hObject, eventdata, handles)
% hObject    handle to EF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global f_span
global centre_f
global start_f
global end_f
start_f = str2double(get(hObject,'String'));
end_f = str2double(get(handles.EF,'String'));
%ef = get(handles.EF,'Value')
f_span = end_f-start_f;
centre_f = start_f + f_span/2;
set(handles.CF,'String',num2str(centre_f));
set(handles.FS,'String',num2str(f_span));

% Hints: get(hObject,'String') returns contents of EF as text
%        str2double(get(hObject,'String')) returns contents of EF as a double


% --- Executes during object creation, after setting all properties.
function SF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function EF_Callback(hObject, eventdata, handles)
% hObject    handle to EF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global f_span
global centre_f
global start_f
global end_f
start_f = str2double(get(handles.SF,'String'));
end_f = str2double(get(handles.EF,'String'));
%ef = get(handles.EF,'Value')
f_span = end_f-start_f;
centre_f = start_f + f_span/2;
set(handles.CF,'String',num2str(centre_f));
set(handles.FS,'String',num2str(f_span));
% Hints: get(hObject,'String') returns contents of EF as text
%        str2double(get(hObject,'String')) returns contents of EF as a double


% --- Executes during object creation, after setting all properties.
function EF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function CF_Callback(hObject, eventdata, handles)
% hObject    handle to EF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global f_span
global centre_f
global start_f
global end_f
f_span = str2double(get(handles.FS,'String'));
centre_f = str2double(get(handles.CF,'String'));
if centre_f<24
    msgbox('Please enter  center frequency in the range of 24 MHz to 1766 MHz ','Error','error','modal');
else if centre_f>1766
        msgbox('Please enter  center frequency in the range of 24 MHz to 1766 MHz ','Error','error','modal');
    else
%ef = get(handles.EF,'Value')
        start_f = centre_f - f_span/2;
        end_f = centre_f + f_span/2;
        set(handles.SF,'String',num2str(start_f));
        set(handles.EF,'String',num2str(end_f));
    end
end
% Hints: get(hObject,'String') returns contents of EF as text
%        str2double(get(hObject,'String')) returns contents of EF as a double


% --- Executes during object creation, after setting all properties.
function CF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end







% --- Executes on button press in Auto_Gain.
function Auto_Gain_Callback(hObject, eventdata, handles)
% hObject    handle to Auto_Gain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(hObject,'Value')==1
    set(handles.Gain,'Enable','Off');
else set(handles.Gain,'Enable','On');
end
% Hint: get(hObject,'Value') returns toggle state of Auto_Gain


% --- Executes on selection change in Sample_Frame.
function Sample_Frame_Callback(hObject, eventdata, handles)
% hObject    handle to Sample_Frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns Sample_Frame contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Sample_Frame


% --- Executes during object creation, after setting all properties.
function Sample_Frame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Sample_Frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function FS_Callback(hObject, eventdata, handles)
% hObject    handle to FS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global f_span
global centre_f
global start_f
global end_f
f_span = str2double(get(handles.FS,'String'));
centre_f = str2double(get(handles.CF,'String'));
if f_span <=0.225
    msgbox('Please enter a value larger than 0.225 MHz','Error','error','modal');
else 
     start_f = centre_f - f_span/2;
     end_f = centre_f + f_span/2;
     set(handles.SF,'String',num2str(start_f));
     set(handles.EF,'String',num2str(end_f));
end
if ((f_span > 0.225) && (f_span <= 0.3))||((f_span > 0.9)&&(f_span<=2.56))
    start_f = centre_f - f_span/2;
     end_f = centre_f + f_span/2;
     set(handles.SF,'String',num2str(start_f));
     set(handles.EF,'String',num2str(end_f));
     set(handles.FBS,'String',num2str(f_span*1000));
     set(handles.FBS,'Enable','Off');
     set(handles.BO,'String','0');
     set(handles.BO,'Enable','Off');
else
    set(handles.FBS,'Enable','On');
    set(handles.BO,'String','35');
    set(handles.BO,'Enable','On');
    
end

% Hints: get(hObject,'String') returns contents of FS as text
%        str2double(get(hObject,'String')) returns contents of FS as a double


% --- Executes on button press in Hardware_Report.
function Hardware_Report_Callback(hObject, eventdata, handles)
% hObject    handle to Hardware_Report (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.SF,'enable','Off');
set(handles.EF,'enable','Off');
set(handles.CF,'enable','Off');
set(handles.FS,'enable','Off');
set(handles.display,'enable','Off');
set(handles.Auto_Scale,'enable','Off');
set(handles.Smooth,'enable','Off');
set(handles.Y_MAX,'enable','Off');
set(handles.Y_MIN,'enable','Off');
set(handles.PH,'enable','Off');
set(handles.OPH,'enable','Off');
set(handles.AVG,'enable','Off');
set(handles.OPH,'enable','Off');
set(handles.Sample_Frame,'enable','Off');
set(handles.Windowing,'enable','Off');
set(handles.Auto_Gain,'enable','Off');
set(handles.Gain,'enable','Off');
set(handles.Hardware_Report,'enable','Off');
set(handles.Start,'enable','Off');
set(handles.Stop,'enable','Off');
set(handles.FBS,'Enable','Off');
set(handles.BO,'Enable','Off');
set(handles.Path,'enable','Off');
set(handles.BaseName,'enable','Off');
set(handles.SweepPerFile,'enable','Off');
set(handles.FileType,'enable','Off');
set(handles.Browse,'enable','Off');
set(handles.Save_data,'enable','Off');
set(handles.Remove_DC,'enable','Off');
set(handles.SweepsPerPlot,'Enable','Off');
set(handles.PSD,'enable','Off');
set(handles.TwoD,'enable','Off');
set(handles.ThreeD,'enable','Off');
set(handles.Credit,'enable','Off');
drawnow

info=sdrinfo();
a = ['RadioName: ',info.RadioName];
b = ['RadioAddress: ',info.RadioAddress];
c = ['RadioIsOpen: ',num2str(info.RadioIsOpen)];
d = ['TunerName: ',info.TunerName];
e = ['Manufacturer: ',info.Manufacturer];
f = ['Product: ',info.Product];
s = num2str(length(info.GainValues));
g = ['GainValues: ','[',s,'x1 double]'];
h = ['RTLCrystalFrequency: ',num2str(info.RTLCrystalFrequency)];
i = ['TunerCrystalFrequency: ', num2str(info.TunerCrystalFrequency)];
j = ['SamplingMode: ',info.SamplingMode];
k = ['OffsetTuning: ',info.OffsetTuning];
waitfor(msgbox({a;b;c;d;e;f;g;h;i;j;k},'Hardware Information', 'modal'));
set(handles.SF,'enable','On');
set(handles.EF,'enable','On');
set(handles.CF,'enable','On');
set(handles.FS,'enable','On');
set(handles.display,'enable','On');
set(handles.Auto_Scale,'enable','On');
set(handles.Smooth,'enable','On');
set(handles.Y_MAX,'enable','On');
set(handles.Y_MIN,'enable','On');
set(handles.Sample_Frame,'enable','On');
set(handles.Windowing,'enable','On');
set(handles.Auto_Gain,'enable','On');
set(handles.Hardware_Report,'enable','On');
set(handles.Start,'enable','On');
set(handles.Stop,'enable','On');
set(handles.Credit,'enable','On');
if str2double(get(handles.FS,'string'))*1e3~=str2double(get(handles.FBS,'string'))
   set(handles.FBS,'Enable','On');
   set(handles.BO,'Enable','On');
end
set(handles.Remove_DC,'enable','On');
if get(handles.Auto_Gain,'Value')==0
   set(handles.Gain,'enable','On');
end
set(handles.Save_data,'enable','On');
if get(handles.Save_data,'Value')==1
    set(handles.Path,'enable','On');
    set(handles.BaseName,'enable','On');
    set(handles.SweepPerFile,'enable','On');
    set(handles.FileType,'enable','On');
    set(handles.Browse,'enable','On');
end
set(handles.PSD,'enable','On');
set(handles.TwoD,'enable','On');
set(handles.ThreeD,'enable','On');
if get(handles.PSD,'Value')==1
    
    set(handles.PH,'Enable','On');
    set(handles.OPH,'Enable','On');
    set(handles.AVG,'Enable','On');
else set(handles.SweepsPerPlot,'Enable','On')
end

function Gain_Callback(hObject, eventdata, handles)
% hObject    handle to Gain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Gain as text
%        str2double(get(hObject,'String')) returns contents of Gain as a double


% --- Executes during object creation, after setting all properties.
function Gain_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Gain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Smooth.
function Smooth_Callback(hObject, eventdata, handles)
% hObject    handle to Smooth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global length_data
global data_max
global data_min
global data_avg
data_max = linspace(0,0,length_data); 
data_min = linspace(0,0,length_data);
data_avg = linspace(0,0,length_data);
% Hint: get(hObject,'Value') returns toggle state of Smooth



function FBS_Callback(hObject, eventdata, handles)
% hObject    handle to FBS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global f_span
fbs = str2double(get(hObject,'String'));


    if f_span <= 0.3
        if fbs > f_span*1000
           msgbox(['Please enter the recommend value ', num2str(f_span*1000),' KHz or enter the value in the range of (225,',num2str(f_span*1000),'] KHz'],'Error','error','modal');
        else if fbs <= 225
                msgbox(['Please enter the value in the range of (225,',num2str(f_span*1000),'] KHz'],'Error','error','modal');
             end
        end
    end
    if (f_span > 0.3) && (f_span <= 0.9)
        if fbs > f_span*1000
           msgbox('Please enter the value in the range of (225,300]','Error','error','modal');
        else if fbs > 300
                msgbox('Please enter the value in the range of (225,300]','Error','error','modal');
             else if fbs <= 225
                      msgbox('Please enter the value in the range of (225,300]','Error','error','modal');
                 end
             end
        end 
    end
    if (f_span > 0.9) && (f_span <= 2.56)
        if fbs>f_span*1000
            msgbox(['Please enter the value in the range of (225,300] and (900,',num2str(f_span*1000),'] KHz'],'Error','error','modal');
        else if fbs <=900
                if fbs > 300
                    msgbox(['Please enter the value in the range of (225,300] and (900,',num2str(f_span*1000),'] KHz'],'Error','error','modal');
                else if fbs <= 225
                        msgbox(['Please enter the value in the range of (225,300] and (900,',num2str(f_span*1000),'] KHz'],'Error','error','modal');
                     end
                end
             end
        end
    end
    if f_span > 2.56
        if fbs>f_span*1000
            msgbox('Please enter the value in the range of (225,300] and (900,2560] KHz','Error','error','modal');
        else if fbs > 2560
                msgbox('Please enter the value in the range of (225,300] and (900,2560] KHz','Error','error','modal');
            
             else if fbs <=900
                     if fbs > 300
                         msgbox(['Please enter the value in the range of (225,300] and (900,',num2str(f_span*1000),'] KHz'],'Error','error','modal');
                     else if fbs <= 225
                             msgbox(['Please enter the value in the range of (225,300] and (900,',num2str(f_span*1000),'] KHz'],'Error','error','modal');
                         end
                     end
                 end
            end
        end
   end


% Hints: get(hObject,'String') returns contents of FBS as text
%        str2double(get(hObject,'String')) returns contents of FBS as a double


% --- Executes during object creation, after setting all properties.
function FBS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FBS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function BO_Callback(hObject, eventdata, handles)
% hObject    handle to BO (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BO as text
%        str2double(get(hObject,'String')) returns contents of BO as a double


% --- Executes during object creation, after setting all properties.
function BO_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BO (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global h
release(h);

% Hint: delete(hObject) closes the figure
delete(hObject);


% --- Executes on button press in Remove_DC.
function Remove_DC_Callback(hObject, eventdata, handles)
% hObject    handle to Remove_DC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Remove_DC


% --- Executes on button press in Save_data.
function Save_data_Callback(hObject, eventdata, handles)
% hObject    handle to Save_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global basename
global sweepperfile
global filetype
global filetype_num
global path
path = get(handles.Path,'string');
basename = get(handles.BaseName,'string');
sweepperfile = str2double(get(handles.SweepPerFile,'string'));
contents = cellstr(get(handles.FileType,'String'));
filetype_num = get(handles.FileType,'Value');
filetype = contents{filetype_num};
if get(hObject,'Value')==1
    set(handles.Path,'enable','On');
    set(handles.BaseName,'enable','On');
    set(handles.SweepPerFile,'enable','On');
    set(handles.FileType,'enable','On');
    set(handles.Browse,'enable','On');
else
    set(handles.Path,'enable','Off');
    set(handles.BaseName,'enable','Off');
    set(handles.SweepPerFile,'enable','Off');
    set(handles.FileType,'enable','Off');
    set(handles.Browse,'enable','Off');
end
% Hint: get(hObject,'Value') returns toggle state of Save_data



function Path_Callback(hObject, eventdata, handles)
% hObject    handle to Path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Path as text
%        str2double(get(hObject,'String')) returns contents of Path as a double


% --- Executes during object creation, after setting all properties.
function Path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function BaseName_Callback(hObject, eventdata, handles)
% hObject    handle to BaseName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BaseName as text
%        str2double(get(hObject,'String')) returns contents of BaseName as a double


% --- Executes during object creation, after setting all properties.
function BaseName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BaseName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SweepPerFile_Callback(hObject, eventdata, handles)
% hObject    handle to SweepPerFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SweepPerFile as text
%        str2double(get(hObject,'String')) returns contents of SweepPerFile as a double


% --- Executes during object creation, after setting all properties.
function SweepPerFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SweepPerFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in FileType.
function FileType_Callback(hObject, eventdata, handles)
% hObject    handle to FileType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns FileType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from FileType


% --- Executes during object creation, after setting all properties.
function FileType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FileType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Browse.
function Browse_Callback(hObject, eventdata, handles)
% hObject    handle to Browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global path
path = uigetdir;
set(handles.Path,'string',path);



function SweepsPerPlot_Callback(hObject, eventdata, handles)
% hObject    handle to SweepsPerPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global mean_time
global psd
global spectrogram_x_axis
global spectrogram_y_axis
global SweepsPerPlot
global SampleRate
global start_f
global end_f
global length_data
global f_span
global off
SweepsPerPlot=str2double(get(hObject,'String'));
if off == 0
    psd=zeros(length_data,SweepsPerPlot);
    for n=1:SweepsPerPlot

        if f_span*1e6 ~= SampleRate
            spectrogram_x_axis(n,:)=linspace(start_f, (start_f*1e6 + ceil((end_f - start_f)*1e6/SampleRate)*SampleRate)/1e6, length_data);
            spectrogram_y_axis(n,:) = linspace(n*mean_time,n*mean_time,length_data);
        else
            spectrogram_x_axis(n,:) = linspace(start_f,end_f,length_data);
            spectrogram_y_axis(n,:) = linspace(n*mean_time,n*mean_time,length_data);

        end

    end
end

% Hints: get(hObject,'String') returns contents of SweepsPerPlot as text
%        str2double(get(hObject,'String')) returns contents of SweepsPerPlot as a double


% --- Executes during object creation, after setting all properties.
function SweepsPerPlot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SweepsPerPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PSD.
function PSD_Callback(hObject, eventdata, handles)
% hObject    handle to PSD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global reduction
global plot_type 
plot_type = 1;
reduction = 1;
set(handles.SweepsPerPlot,'Enable','Off');
set(handles.PH,'Enable','On');
set(handles.OPH,'Enable','On');
set(handles.AVG,'Enable','On');
set(handles.Reduc,'Enable','Off');
set(handles.Reduc,'Value',1);
set(handles.Reduct_text,'string',['1:',num2str(reduction)]);
% Hint: get(hObject,'Value') returns toggle state of PSD


% --- Executes on button press in TwoD.
function TwoD_Callback(hObject, eventdata, handles)
% hObject    handle to TwoD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global waterfall_loop
global psd
global length_data
global SweepsPerPlot

global reduction
global plot_type 

plot_type = 2;
reduction = 1;
waterfall_loop = 1;
psd=zeros(length_data,SweepsPerPlot);

set(handles.SweepsPerPlot,'Enable','On');

set(handles.PH,'Enable','Off');
set(handles.OPH,'Enable','Off');
set(handles.AVG,'Enable','Off');
set(handles.Reduc,'Enable','Off');
set(handles.Reduc,'Value',1);
set(handles.PH,'value',0);
set(handles.OPH,'value',0);
set(handles.AVG,'value',0);
set(handles.Reduct_text,'string',['1:',num2str(reduction)]);
% Hint: get(hObject,'Value') returns toggle state of TwoD


% --- Executes on button press in ThreeD.
function ThreeD_Callback(hObject, eventdata, handles)
% hObject    handle to ThreeD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global waterfall_loop
global psd
global length_data
global SweepsPerPlot
global plot_type 
plot_type = 4;
if get(handles.display,'value')==2
    psd=(-100).*ones(length_data,SweepsPerPlot);
else
    psd=zeros(length_data,SweepsPerPlot);
end

waterfall_loop = 1;



set(handles.SweepsPerPlot,'Enable','On');

set(handles.PH,'Enable','Off');
set(handles.OPH,'Enable','Off');
set(handles.AVG,'Enable','Off');
set(handles.PH,'value',0);
set(handles.OPH,'value',0);
set(handles.AVG,'value',0);
set(handles.Reduc,'Enable','On');

% Hint: get(hObject,'Value') returns toggle state of ThreeD


% --- Executes on button press in Credit.
function Credit_Callback(hObject, eventdata, handles)
% hObject    handle to Credit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.SF,'enable','Off');
set(handles.EF,'enable','Off');
set(handles.CF,'enable','Off');
set(handles.FS,'enable','Off');
set(handles.FBS,'Enable','Off');
set(handles.BO,'Enable','Off');
set(handles.PSD,'enable','Off');
set(handles.TwoD,'enable','Off');
set(handles.ThreeD,'enable','Off');
set(handles.SweepsPerPlot,'Enable','Off');
set(handles.Remove_DC,'enable','Off');
set(handles.Smooth,'enable','Off');
set(handles.Auto_Gain,'enable','Off');
set(handles.Gain,'enable','Off');
set(handles.PH,'enable','Off');
set(handles.AVG,'enable','Off');
set(handles.OPH,'enable','Off');
set(handles.display,'enable','Off');
set(handles.Y_MIN,'enable','Off');
set(handles.Y_MAX,'enable','Off');
set(handles.Auto_Scale,'enable','Off');
set(handles.Sample_Frame,'enable','Off');
set(handles.Windowing,'enable','Off');
set(handles.Save_data,'enable','Off');
set(handles.BaseName,'enable','Off');
set(handles.SweepPerFile,'enable','Off');
set(handles.FileType,'enable','Off');
set(handles.Path,'enable','Off');
set(handles.Browse,'enable','Off');
set(handles.Start,'enable','Off');
set(handles.Stop,'enable','Off');
set(handles.Hardware_Report,'enable','Off');
set(handles.Credit,'enable','Off');
drawnow

waitfor(msgbox({'Developer: Yibo Zhang';...
                'Department of Electrical Engineering & Electronics';...
                'University of Liverpool';...
                'Liverpool L69 3GJ';...
                'E-mail: sgyzha54@liv.ac.uk';...
                },'Information', 'modal'));
set(handles.SF,'enable','On');
set(handles.EF,'enable','On');
set(handles.CF,'enable','On');
set(handles.FS,'enable','On');
set(handles.display,'enable','On');
set(handles.Auto_Scale,'enable','On');
set(handles.Smooth,'enable','On');
set(handles.Y_MAX,'enable','On');
set(handles.Y_MIN,'enable','On');
set(handles.Sample_Frame,'enable','On');
set(handles.Windowing,'enable','On');
set(handles.Auto_Gain,'enable','On');
set(handles.Hardware_Report,'enable','On');
set(handles.Credit,'enable','On');
set(handles.Start,'enable','On');
set(handles.Stop,'enable','On');
if str2double(get(handles.FS,'string'))*1e3~=str2double(get(handles.FBS,'string'))
   set(handles.FBS,'Enable','On');
   set(handles.BO,'Enable','On');
end
set(handles.Remove_DC,'enable','On');
if get(handles.Auto_Gain,'Value')==0
   set(handles.Gain,'enable','On');
end
set(handles.Save_data,'enable','On');
if get(handles.Save_data,'Value')==1
    set(handles.Path,'enable','On');
    set(handles.BaseName,'enable','On');
    set(handles.SweepPerFile,'enable','On');
    set(handles.FileType,'enable','On');
    set(handles.Browse,'enable','On');
end
set(handles.PSD,'enable','On');
set(handles.TwoD,'enable','On');
set(handles.ThreeD,'enable','On');
if get(handles.PSD,'Value')==1
    
    set(handles.PH,'Enable','On');
    set(handles.OPH,'Enable','On');
    set(handles.AVG,'Enable','On');
else set(handles.SweepsPerPlot,'Enable','On')
end


% --- Executes on slider movement.
function Reduc_Callback(hObject, eventdata, handles)
% hObject    handle to Reduc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global reduction
reduction = ceil(get(hObject,'Value'));
set(handles.Reduct_text,'string',['1:',num2str(reduction)]);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function Reduc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Reduc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes during object creation, after setting all properties.
function Reduct_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Reduct_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in test.
function test_Callback(hObject, eventdata, handles)
% hObject    handle to test (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
f = figure;
copyobj(handles.axes1,f);
