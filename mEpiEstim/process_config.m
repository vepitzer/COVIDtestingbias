function config = process_config(config) 

% THIS FUNCTION IS NOT NEEDED. ALL CHECKS ARE DONE IN CLASS

    if isempty(config.mean_prior)
        config.mean_prior = 5;
    end

    if isempty(config.std_prior)
        config.std_prior = 5;
    end

%     if config.mean_prior <= 0
%         error("config$mean_prior must be >0.")
%     end
%     if config.std_prior <= 0
%         error("config$std_prior must be >0.")
%     end

    if isempty(config.cv_posterior)
        config.cv_posterior = 0.3;
    end
    
    if isempty(config.mcmc_control)
        config.mcmc_control = make_mcmc_control();
    end
    
end


