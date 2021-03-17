clc;
close all;
clear variables;

load('IS.mat');


baselineValid = datatimes > eventTimes(1) - seconds(120) & datatimes < eventTimes(1) - seconds(avSampleTime);
baseline = data(:,baselineValid);
baseline_bg = bg_data(:,baselineValid);
baseline_fg = fg_data(:,baselineValid);
baseline_av = mean(baseline,2);
baseline_bg_av = mean(baseline_bg,2);
baseline_fg_av = mean(baseline_fg,2);

baseline_tot_av = mean(nansum(baseline,1));
baseline_tot_bg_av = mean(nansum(baseline_bg,1));
baseline_tot_fg_av = mean(nansum(baseline_fg,1));

breathingValid_nodev = datatimes > eventTimes(1) & datatimes < eventTimes(2)+seconds(avSampleTime);
breathing_nodev = data(:,breathingValid_nodev);
breathing_bg_nodev = bg_data(:,breathingValid_nodev);
breathing_fg_nodev = fg_data(:,breathingValid_nodev);
breathing_nodev_times = datatimes(breathingValid_nodev);

breathingValid_dev = datatimes > eventTimes(7) & datatimes < eventTimes(8)+seconds(avSampleTime);
breathing_dev = data(:,breathingValid_dev);
breathing_bg_dev = bg_data(:,breathingValid_dev);
breathing_fg_dev = fg_data(:,breathingValid_dev);
breathing_dev_times = datatimes(breathingValid_dev);

speakingValid_nodev = datatimes > eventTimes(3) & datatimes < eventTimes(4)+seconds(avSampleTime);
speaking_nodev = data(:,speakingValid_nodev);
speaking_bg_nodev = bg_data(:,speakingValid_nodev);
speaking_fg_nodev = fg_data(:,speakingValid_nodev);
speaking_nodev_times = datatimes(speakingValid_nodev);

speakingValid_dev = datatimes > eventTimes(9) & datatimes < eventTimes(10)+seconds(avSampleTime);
speaking_dev = data(:,speakingValid_dev);
speaking_bg_dev = bg_data(:,speakingValid_dev);
speaking_fg_dev = fg_data(:,speakingValid_dev);
speaking_dev_times = datatimes(speakingValid_dev);

coughingValid_nodev = datatimes > eventTimes(5) & datatimes < eventTimes(6)+seconds(avSampleTime);
coughing_nodev = data(:,coughingValid_nodev);
coughing_bg_nodev = bg_data(:,coughingValid_nodev);
coughing_fg_nodev = fg_data(:,coughingValid_nodev);
coughing_nodev_times = datatimes(coughingValid_nodev);

coughingValid_dev = datatimes > eventTimes(11) & datatimes < eventTimes(12)+seconds(avSampleTime);
coughing_dev = data(:,coughingValid_dev);
coughing_bg_dev = bg_data(:,coughingValid_dev);
coughing_fg_dev = fg_data(:,coughingValid_dev);
coughing_dev_times = datatimes(coughingValid_dev);

nSizes = size(data,1);
tColor = lines(nSizes);

for m=1:2
    
    if m==1
        time = breathing_nodev_times - breathing_nodev_times(1);
        data_temp = breathing_nodev;
        bg_temp = breathing_bg_nodev;
        fg_temp = breathing_fg_nodev;
    elseif m==2
        time = breathing_dev_times - breathing_dev_times(1);
        data_temp = breathing_dev;
        bg_temp = breathing_bg_dev;
        fg_temp = breathing_fg_dev;
    end
    
%     if m==1
%         time = speaking_nodev_times - speaking_nodev_times(1);
%         data_temp = speaking_nodev;
%         bg_temp = speaking_bg_nodev;
%         fg_temp = speaking_fg_nodev;
%     elseif m==2
%         time = speaking_dev_times - speaking_dev_times(1);
%         data_temp = speaking_dev;
%         bg_temp = speaking_bg_dev;
%         fg_temp = speaking_fg_dev;
%     end
%     
%         
%     if m==1
%         time = coughing_nodev_times - coughing_nodev_times(1);
%         data_temp = coughing_nodev;
%         bg_temp = coughing_bg_nodev;
%         fg_temp = coughing_fg_nodev;
%     elseif m==2
%         time = coughing_dev_times - coughing_dev_times(1);
%         data_temp = coughing_dev;
%         bg_temp = coughing_bg_dev;
%         fg_temp = coughing_fg_dev;
%     end

    for k=1:nSizes+1

        if (k <= nSizes)
            subplot(nSizes+1,3,3*(k-1)+1);
            if m==1
                plot(time, data_temp(k,:),'--','Color',tColor(k,:));
                yline(baseline_av(k),'k:', 'LineWidth', 2);
            else
                plot(time, data_temp(k,:),'Color',tColor(k,:));
            end

            title(['Diameter: ', num2str(diameters(k)), '\mum']);
            ylabel('#/m^3');
            xlabel('time');
            xline(seconds(0),'k:');
            hold on;
            ylim_temp = ylim;
            ylim([0, max(ylim_temp)]);

            subplot(nSizes+1,3,3*(k-1)+2);
            if m==1
                plot(time, bg_temp(k,:),'--','Color',tColor(k,:));
                yline(baseline_bg_av(k),'k:');
            else
                plot(time, bg_temp(k,:),'Color',tColor(k,:));
            end

            title(['Diameter: ', num2str(diameters(k)), '\mum']);
            ylabel('#/m^3');
            xlabel('time');
            xline(seconds(0),'k:');
            hold on;
            ylim_temp = ylim;
            ylim([0, max(ylim_temp)]);

            subplot(nSizes+1,3,3*(k-1)+3);
            if m==1
                plot(time, fg_temp(k,:),'--','Color',tColor(k,:));
                %legend('without device');
                yline(baseline_bg_av(k),'k:');
            else
                plot(time, fg_temp(k,:),'Color',tColor(k,:));
            end

            title(['Diameter: ', num2str(diameters(k)), '\mum']);
            ylabel('#/m^3');
            xlabel('time');
            xline(seconds(0),'k:');
            hold on;
            ylim_temp = ylim;
            ylim([0, max(ylim_temp)]);
            
            %pause(0.1);
        else
            if k== nSizes+1
                subplot(nSizes+1,3,3*(k-1)+1);
                if m==1
                    plot(time, nansum(data_temp,1),'k--');
                    yline(baseline_tot_av,'k:', 'LineWidth', 2);
                else
                    plot(time,nansum(data_temp,1),'k');
                end

                title(['Total #']);
                ylabel('#/m^3');
                xlabel('time');
                xline(seconds(0),'k:');
                hold on;
                ylim_temp = ylim;
                ylim([0, max(ylim_temp)]);

                subplot(nSizes+1,3,3*(k-1)+2);
                if m==1
                    plot(time, nansum(bg_temp,1),'k--');
                    yline(baseline_tot_bg_av,'k:', 'LineWidth', 2);
                else
                    plot(time,nansum(bg_temp,1),'k');
                end

                title(['Total #']);
                ylabel('#/m^3');
                xlabel('time');
                xline(seconds(0),'k:');
                hold on;
                ylim_temp = ylim;
                ylim([0, max(ylim_temp)]);

                subplot(nSizes+1,3,3*(k-1)+3);
                if m==1
                    plot(time, nansum(fg_temp,1),'k--');
                    yline(baseline_tot_fg_av,'k:', 'LineWidth', 2);
                else
                    plot(time,nansum(fg_temp,1),'k');
                end

                title(['Total #']);
                ylabel('#/m^3');
                xlabel('time');
                xline(seconds(0),'k:');
                hold on;
                ylim_temp = ylim;
                ylim([0, max(ylim_temp)]);
            elseif k == nSizes+2

            end

        end
    end
end


speakingValid_nodev = datatimes > eventtimes(3) && datatimes < eventtimes(4);
