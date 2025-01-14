function img=hnd_read(filename);
% function img=read_hndFile(filename);

if(nargin<1)
    error([mfilename,' requires 1 parameter: filename.']);
    return
end

if (~exist(filename))
  error([mfilename, ': could not find ', filename]);
  return
end

fid = fopen(filename, 'r');
if fid==-1
  error(['Error while opening file ', filename, ' for reading.']);
  return
end
total=0;
bitHist=[0,0,0];
%HNA Header information read here.
img.sFileType=char(fread(fid, 32, 'char'))';        %the file type; defined values are:
total=total+32;                                      %VARIAN_VA_INTERNAL_1.0 or VARIAN_VA_INTERNAL_DYN_COMP_1.0
img.lFilelength=fread(fid, 1, 'uint32');            %file length in bytes
total=total+4;
img.sChecksumSpec=char(fread(fid, 4 , 'char'))';    %checksum specifier; define values are: NONE, ZIP
total=total+4;
img.lCheckSum=fread(fid, 1, 'int32');               %file checksum, always 0
total=total+4;
img.sCreationDate=char(fread(fid, 8, 'char'))';     %the creation date of the image; format: YYYYMMDD
total=total+8;
img.sCreationTime=char(fread(fid, 8, 'char'))';     %the creation time of the image; format: HH:MM:SS
total=total+8;
img.sPatientId=char(fread(fid, 16, 'char'))';       %the patient ID
total=total+16;
img.lPatientSerNr=fread(fid, 1, 'int32');           %Patient Serial Number
total=total+4;
img.sSeriesId=char(fread(fid, 16, 'char'))';        %Series ID
total=total+16;
img.lSeriesSerNr=fread(fid, 1, 'int32');            %Series Serial Number
total=total+4;
img.sSliceId=char(fread(fid, 16, 'char'))';         %Slice ID
total=total+16;
img.lSliceSerNr=fread(fid, 1, 'int32');             %Series Serial Number
total=total+4;
img.lXSize=fread(fid, 1, 'int32');                  %number of columns of the slice (width)
total=total+4;
img.lYSize=fread(fid, 1, 'int32');                  %number of rows of the slice (height)                           
total=total+4;
img.dSliceZPos=fread(fid, 1, 'float64');            %the z position of the slice
total=total+8;
img.sModality=char(fread(fid, 16, 'char'))';        %the modality of the slice; defined values are:
total=total+16;                                     %'CR', 'CT', 'NM', 'MR', 'US'
img.lWindow=fread(fid, 1, 'int32');                 %window/level values
total=total+4;
img.lLevel=fread(fid, 1, 'int32');
total=total+4;
%End reading HNA header info
%Start reading HNC Header info
img.lPixelOffset=fread(fid, 1, 'int32');
total=total+4;
img.szImageType=char(fread(fid, 4, 'char'))';
total=total+4;
MaxMachineInfo=40;

img.MachineAxisGantryRtn=fread(fid, 1, 'float64');      %Gantry rotation angle, clockwise viewed from isocenter. 0° is up.
total=total+8;
img.MachineAxisSAD=fread(fid, 1, 'float64');            %Source Axis Distance, radiation source to isocenter, 100.0 for Clinac, Acuity 
total=total+8;
img.MachineAxisSFD=fread(fid, 1, 'float64');            %Source Film Distance, AcqII specific
total=total+8;
img.MachineAxisCollX1=fread(fid, 1, 'float64');         %Collimator right X position viewed from radiation source. Values increase to right.
total=total+8;
img.MachineAxisCollX2=fread(fid, 1, 'float64');         %Collimator left X position. Values increase to right.
total=total+8;
img.MachineAxisCollY1=fread(fid, 1, 'float64');         %Collimator top (farther from gantry) Y position. Values increase when moving away from gantry.
total=total+8;
img.MachineAxisCollY2=fread(fid, 1, 'float64');         %Collimator bottom Y position. Values increase when moving away from gantry.
total=total+8;
img.MachineAxisCollRtn=fread(fid, 1, 'float64');        %Collimator rotation angle clockwise viewed from isocenter. 0° is toward gantry.
total=total+8;
img.MachineAxisFieldX=fread(fid, 1, 'float64');         %Field
total=total+8;
img.MachineAxisFieldY=fread(fid, 1, 'float64');
total=total+8;
img.MachineAxisBladeX1=fread(fid, 1, 'float64');        %Blades
total=total+8;
img.MachineAxisBladeX2=fread(fid, 1, 'float64');
total=total+8;
img.MachineAxisBladeY1=fread(fid, 1, 'float64');
total=total+8;
img.MachineAxisBladeY2=fread(fid, 1, 'float64');
total=total+8;
img.MachineAxisIDUPosLng=fread(fid, 1, 'float64');      %Image Detection Unit support arm longitudinal position in respect to the isocenter.
total=total+8;
img.MachineAxisIDUPosLat=fread(fid, 1, 'float64');      %IDU lateral position in respect to the isocenter.
total=total+8;
img.MachineAxisIDUPosVrt=fread(fid, 1, 'float64');      %IDU vertical position in respect to the isocenter. SID = SAD - IDUPosVrt
total=total+8;
img.MachineAxisIDUPosRtn=fread(fid, 1, 'float64');      %IDU rotation angle. Always 0.0 for II, IAS2, IAS3.
total=total+8;
img.MachineAxisPatientSupportAngle=fread(fid, 1, 'float64');    %Couch patient support angle = couch isocentric rotation.
total=total+8;                                                  %Values increase counter-clockwise when viewed from above.
img.MachineAxisTableTopEccentricAngle=fread(fid, 1, 'float64'); %Couch table top eccentric angle = couch eccentric rotation.
total=total+8;
img.MachineAxisCouchVrt=fread(fid, 1, 'float64');       %Couch vertical position. Values increase when moving up, 0 is at isocenter.
total=total+8;
img.MachineAxisCouchLng=fread(fid, 1, 'float64');       %Couch longitudinal position. Value increases when moving toward gantry.
total=total+8;
img.MachineAxisCouchLat=fread(fid, 1, 'float64');       %Couch lateral position. Values increase to left when viewed from gantry.
total=total+8;
%Opposite to the machine axis positions, the detector resolutions are in
%millimeters!!!
img.MachineIDUResolutionX=fread(fid, 1, 'float64');     %Horizontal distance [mm] between two pixels on detector plane
total=total+8;
img.MachineIDUResolutionY=fread(fid, 1, 'float64');     %Vertical distance [mm] between two pixels on detector plane
total=total+8;
img.MachineImageResolutionX=fread(fid, 1, 'float64');   %Horizontal distance [mm] between two pixels in respect to the iso centre
total=total+8;
img.MachineImageResolutionY=fread(fid, 1, 'float64');   %Vertical distance [mm] between two pixels in respect to the iso centre
%ImageResolutionX = IDUResolutionX * SAD / (SAD - IDUPosVrt)
total=total+8;
img.MachineEnergy=fread(fid, 1, 'float64');             %Clinac energy [MV]
total=total+8;
img.MachineDoseRate=fread(fid, 1, 'float64');           %Clinac doserate [MU/Min]
total=total+8;
img.MachineXRayKV=fread(fid, 1, 'float64');             %Diagnostic imaging energy [KV]
total=total+8;
img.MachineXRayMA=fread(fid, 1, 'float64');             % Diagnostic imaging doserate [mA]
total=total+8;

%MetersetExposure is the image acquisition duration for one image [min], 
% resp. treatment machine meterset duration over which image has been acquired.
% The acquisition duration depends on the frame rate and the number of averaged frames.
% Simulator last image hold: = beamOffTime - beamOnTime
% Simulator Cine: = always GetNotDefinedDouble
img.MachineMetersetExposure=fread(fid, 1, 'float64');
total=total+8;
img.MachineAcqAdjustment=fread(fid, 1, 'float64');      %Time [min] elapsed since beam on
total=total+8;
% CT and Gating
img.MachineCTProjectionAngle=fread(fid, 1, 'float64');  %[degree]
total=total+8;
img.MachineCTNormChamber=fread(fid, 1, 'float64');
total=total+8;
img.MachineGatingTimeTag=fread(fid, 1, 'float64');
total=total+8;
img.MachineGating4DInfoX=fread(fid, 1, 'float64');
total=total+8;
img.MachineGating4DInfoY=fread(fid, 1, 'float64');
total=total+8;
img.MachineGating4DInfoZ=fread(fid, 1, 'float64');
total=total+8;
img.MachineGating4DInfoTime=fread(fid, 1, 'float64');
total=total+8;
%End reading HNC header info
%Start reading HNC File info
img.szReserved=char(fread(fid, 24, 'char'))';
total=total+24;
%End reading HNC File info
img.dRescaleSlope=fread(fid, 1, 'float64');
total=total+8;
img.dRescaleIntercept=fread(fid, 1, 'float64');
total=total+8;
img.dMinPixelValue=fread(fid, 1, 'float64','b');
total=total+8;
img.dMaxPixelValue=fread(fid, 1, 'float64','l');
total=total+8;

% Varian headder is 1024 so dump the rest of the header that hasn't be
% read.
%disp(total)
total=1024-total;
%disp(total)
%dump=char(fread(fid, total, 'char'))';
fseek(fid,1024,'bof');

%Read the bitTable
bitTablelength=round(img.lXSize*(img.lYSize-1)/4);
bitTable=(fread(fid, bitTablelength, 'uint8'));
pixels=zeros(img.lXSize, img.lYSize);
%Read the uncompressed first row and first pixel of second row.
pixels(:,1)=double(fread(fid, img.lXSize, 'ulong'));
pixels(1,2)=double(fread(fid, 1, 'ulong'));
pixels(1,100)

% totalPixels=img.lXSize*img.lYSize;
% currentPixel=img.lXSize+2; % this is the next pixel to be decompressed


thisByte=1;
thisBit=0;
for atY=2:1:img.lYSize
    if atY==2 
        atX_start=2; 
    else
        atX_start=1;
    end;
    
    for atX=atX_start:1:img.lXSize
        bits=bitget(bitTable(thisByte), thisBit+1:1:thisBit+2);
        thisBit=thisBit+2;
        if(thisBit==8)
            thisBit=0;
            thisByte=thisByte+1;
            if(thisByte>bitTablelength)
                disp('Error in number of bytes in bitTable Reading');
                thisByte
                bitTablelength
                return
            end
        end
        bits=2*bits(2)+bits(1);
        bitHist(bits+1)=bitHist(bits+1)+1;
        switch bits
            case 0
                [diff, count]=fread(fid, 1, 'int8');
            case 1
                [diff, count]=fread(fid, 1, 'int16');
            case 2
                [diff, count]=fread(fid, 1, 'int32');
            otherwise
                disp('Error in bit depth');
                bits
                return
        end
        diff=double(diff);
        if(count==0)
            disp('end of file');
            img.pixels=pixels;
            disp('Row number ');
            disp(round(atY));
            disp('Col number ');
            disp(round(atX));
            disp('Bits');
            disp(bits);
            disp('bitTablelength')
            disp(round(bitTablelength));
            disp('thisByte');
            disp(round(thisByte));
            disp('thisBit');
            disp(round(thisBit));
            bitHist
            whos diff
            return
        end
       
%        pixels(atX, atY)=...
%        diff ... 
%        -pixels(atX-1, atY-1) ... 
%        +pixels(atX-1, atY) ...
%        +pixels(atX, atY-1);
	nIndex=(atY-1)*img.lXSize+atX;
        pixels(atX,atY)= diff - pixels(nIndex-1-img.lXSize)...
			      + pixels(nIndex-1) ...
			      + pixels(nIndex-img.lXSize);       

    end
    %imagesc(pixels)
    %axis image
end
fclose(fid);
%disp('returning');
img.pixels=pixels;
return


