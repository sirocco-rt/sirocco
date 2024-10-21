/* agn.c */
double agn_init(double r, double lum, double alpha, double freqmin, double freqmax, int ioniz_or_extract, double *f);
double emittance_pow(double freqmin, double freqmax, double alpha);
double emittance_bpow(double freqmin, double freqmax, double alpha);
int photo_gen_agn(PhotPtr p, double r, double alpha, double weight, double f1, double f2, int spectype, int istart, int nphot);
/* anisowind.c */
int randwind_thermal_trapping(PhotPtr p, int *nnscat);
/* atomic_extern_init.c */
/* atomicdata.c */
int get_atomic_data(char masterfile[]);
/* atomicdata_init.c */
int init_atomic_data(void);
/* atomicdata_sub.c */
int atomicdata2file(void);
int index_lines(void);
int index_phot_top(void);
int index_inner_cross(void);
void indexx(int n, float arrin[], int indx[]);
int limit_lines(double freqmin, double freqmax);
int check_xsections(void);
double q21(struct lines *line_ptr, double t);
double q12(struct lines *line_ptr, double t);
double a21(struct lines *line_ptr);
double upsilon(int n_coll, double u0);
void skiplines(FILE *fptr, int nskip);
/* bands.c */
int bands_init(int imode, struct xbands *band);
int ion_bands_init(int mode, double freqmin, double freqmax, struct xbands *band);
void check_appropriate_banding(struct xbands *band, int mode);
/* bb.c */
double planck(double t, double freqmin, double freqmax);
double get_rand_pow(double x1, double x2, double alpha);
double get_rand_exp(double alpha_min, double alpha_max);
double integ_planck_d(double alphamin, double alphamax);
int init_integ_planck_d(void);
double planck_d(double alpha, void *params);
double planck_d_2(double alpha, void *params);
double emittance_bb(double freqmin, double freqmax, double t);
double check_freq_max(double freq_max, double temp);
/* bilinear.c */
int bilin(double x[], double x00[], double x01[], double x10[], double x11[], double *f, double *g);
/* brem.c */
double integ_brem(double freq, void *params);
double brem_d(double alpha, void *params);
double get_rand_brem(double freqmin, double freqmax);
/* cdf.c */
int cdf_gen_from_func(CdfPtr cdf, double (*func)(double, void *), double xmin, double xmax, int njumps, double jump[]);
double gen_array_from_func(double (*func)(double, void *), double xmin, double xmax, int pdfsteps);
int cdf_gen_from_array(CdfPtr cdf, double x[], double y[], int n_xy, double xmin, double xmax);
double cdf_get_rand(CdfPtr cdf);
int cdf_limit(CdfPtr cdf, double xmin, double xmax);
double cdf_get_rand_limit(CdfPtr cdf);
int cdf_to_file(CdfPtr cdf, char comment[]);
int cdf_inputs_to_file(double x[], double y[], int n_xy, double xmin, double xmax, char filename[]);
int cdf_check(CdfPtr cdf);
int calc_cdf_gradient(CdfPtr cdf);
int cdf_array_fixup(double *x, double *y, int n_xy);
/* charge_exchange.c */
int compute_ch_ex_coeffs(double T);
double ch_ex_heat(WindPtr one, double t_e);
/* communicate_macro.c */
void broadcast_macro_atom_emissivities(const int n_start, const int n_stop, const int n_cells_rank);
void broadcast_macro_atom_recomb(const int n_start, const int n_stop, const int n_cells_rank);
int broadcast_updated_macro_atom_properties(const int n_start, const int n_stop, const int n_cells_rank);
int broadcast_macro_atom_state_matrix(int n_start, int n_stop, int n_cells_rank);
void reduce_macro_atom_estimators(void);
/* communicate_plasma.c */
void broadcast_plasma_grid(const int n_start, const int n_stop, const int n_cells_rank);
void broadcast_wind_luminosity(const int n_start, const int n_stop, const int n_cells_rank);
void broadcast_wind_cooling(const int n_start, const int n_stop, const int n_cells_rank);
int broadcast_updated_plasma_properties(const int n_start_rank, const int n_stop_rank, const int n_cells_rank);
int reduce_simple_estimators(void);
/* communicate_spectra.c */
int normalize_spectra_across_ranks(void);
/* communicate_wind.c */
void broadcast_wind_grid(const int n_start, const int n_stop, const int n_cells_rank);
/* compton.c */
int compton_scatter(PhotPtr p);
double kappa_comp(PlasmaPtr xplasma, double freq);
double kappa_ind_comp(PlasmaPtr xplasma, double freq);
double total_comp(WindPtr one, double t_e);
double klein_nishina(double nu);
void set_comp_func_values(double rand_cs, double max_cs, double energy_ratio);
int compton_dir(PhotPtr p);
double pdf_thermal(double x, void *params);
int compton_get_thermal_velocity(double t, double *v);
double compton_func(double f, void *params);
double sigma_compton_partial(double f, double x);
double compton_alpha(double nu);
double compton_beta(double nu);
double comp_cool_integrand(double nu, void *params);
double compton_reweight_norm(double nu);
int compton_reweight(PhotPtr p_in, PhotPtr p_out);
/* continuum.c */
double one_continuum(int spectype, double t, double g, double freqmin, double freqmax);
double emittance_continuum(int spectype, double freqmin, double freqmax, double t, double g);
double model_int(double lambda, void *params);
/* cooling.c */
double cooling(PlasmaPtr xplasma, double t);
double xtotal_emission(WindPtr one, double f1, double f2);
double adiabatic_cooling(WindPtr one, double t);
double shock_heating(WindPtr one);
double wind_cooling(void);
/* corona.c */
int get_corona_params(int ndom);
double corona_velocity(int ndom, double x[], double v[]);
double corona_rho(int ndom, double x[]);
/* cv.c */
double wdrad(double m);
double diskrad(double m1, double m2, double period);
double roche2(double q, double a);
double logg(double mass, double rwd);
/* cylind_var.c */
double cylvar_ds_in_cell(int ndom, PhotPtr p);
int cylvar_make_grid(int ndom, WindPtr w);
int cylvar_wind_complete(int ndom, WindPtr w);
int cylvar_cell_volume(WindPtr w);
int cylvar_where_in_grid(int ndom, double x[], int ichoice, double *fx, double *fz);
int cylvar_get_random_location(int n, double x[]);
int cylvar_extend_density(int ndom, WindPtr w);
int cylvar_coord_fraction(int ndom, int ichoice, double x[], int ii[], double frac[], int *nelem);
void cylvar_allocate_domain(int ndom);
/* cylindrical.c */
double cylind_ds_in_cell(int ndom, PhotPtr p);
int cylind_make_grid(int ndom, WindPtr w);
int cylind_wind_complete(int ndom, WindPtr w);
int cylind_cell_volume(WindPtr w);
int cylind_where_in_grid(int ndom, double x[]);
int cylind_get_random_location(int n, double x[]);
int cylind_extend_density(int ndom, WindPtr w);
int cylind_is_cell_in_wind(int n);
/* define_wind.c */
void define_wind(void);
/* density.c */
double get_ion_density(int ndom, double x[], int nion);
/* diag.c */
int get_standard_care_factors(void);
int get_extra_diagnostics(void);
int init_searchlight(void);
int init_extra_diagnostics(void);
int save_photon_stats(WindPtr one, PhotPtr p, double ds, double w_ave);
int save_photons(PhotPtr p, char comment[]);
int track_scatters(PhotPtr p, int nplasma, char *comment);
int Diag(char *format, ...);
/* dielectronic.c */
int compute_dr_coeffs(double temp);
double total_dr(WindPtr one, double t_e);
/* direct_ion.c */
int compute_di_coeffs(double T);
double q_ioniz_dere(int nion, double t_e);
double total_di(WindPtr one, double t_e);
int compute_qrecomb_coeffs(double T);
double q_recomb_dere(struct topbase_phot *cont_ptr, double electron_temperature);
double q_ioniz(struct topbase_phot *cont_ptr, double electron_temperature);
double q_recomb(struct topbase_phot *cont_ptr, double electron_temperature);
/* disk.c */
double teff(double x);
double geff(double x);
double vdisk(double x[], double v[]);
double zdisk(double r);
double ds_to_disk(struct photon *p, int allow_negative, int *hit);
double disk_height(double s, void *params);
double disk_colour_correction(double t);
/* disk_init.c */
double disk_init(double rmin, double rmax, double m, double mdot, double freqmin, double freqmax, int ioniz_or_extract, double *ftot);
int qdisk_init(double rmin, double rmax, double m, double mdot);
int qdisk_reinit(PhotPtr p);
int qdisk_save(char *diskfile, int ichoice);
int read_non_standard_disk_profile(char *tprofile);
/* disk_photon_gen.c */
int photo_gen_disk(PhotPtr p, double weight, double f1, double f2, int spectype, int istart, int nphot);
int disk_photon_summary(char filename[], char mode[]);
/* emission.c */
double wind_luminosity(double f1, double f2, int mode);
double total_emission(PlasmaPtr xplasma, double f1, double f2);
int photo_gen_wind(PhotPtr p, double weight, double freqmin, double freqmax, int photstart, int nphot);
double one_line(PlasmaPtr xplasma, int *nres);
double total_free(PlasmaPtr xplasma, double t_e, double f1, double f2);
double ff(PlasmaPtr xplasma, double t_e, double freq);
double one_ff(PlasmaPtr xplasma, double f1, double f2);
double gaunt_ff(double gsquared);
/* estimators_macro.c */
int bf_estimators_increment(WindPtr one, PhotPtr p, double ds);
int bb_estimators_increment(WindPtr one, PhotPtr p, double tau_sobolev, double dvds, int nn);
int normalise_macro_estimators(PlasmaPtr xplasma);
double total_fb_matoms(PlasmaPtr xplasma, double t_e, double f1, double f2);
double total_bb_cooling(PlasmaPtr xplasma, double t_e);
double macro_bb_heating(PlasmaPtr xplasma, double t_e);
double macro_bf_heating(PlasmaPtr xplasma, double t_e);
int bb_simple_heat(PlasmaPtr xplasma, PhotPtr p, double tau_sobolev, int nn);
int check_stimulated_recomb(PlasmaPtr xplasma);
int get_dilute_estimators(PlasmaPtr xplasma);
double get_gamma(struct topbase_phot *cont_ptr, PlasmaPtr xplasma);
double gamma_integrand(double freq, void *params);
double get_gamma_e(struct topbase_phot *cont_ptr, PlasmaPtr xplasma);
double gamma_e_integrand(double freq, void *params);
double get_alpha_st(struct topbase_phot *cont_ptr, PlasmaPtr xplasma);
double alpha_st_integrand(double freq, void *params);
double get_alpha_st_e(struct topbase_phot *cont_ptr, PlasmaPtr xplasma);
double alpha_st_e_integrand(double freq, void *params);
/* estimators_simple.c */
int update_banded_estimators(PlasmaPtr xplasma, PhotPtr p, double ds, double w_ave, int ndom);
int update_flux_estimators(PlasmaPtr xplasma, PhotPtr phot_mid, double ds_obs, double w_ave, int ndom);
int update_force_estimators(PlasmaPtr xplasma, PhotPtr p, PhotPtr phot_mid, double ds, double w_ave, int ndom, double z, double frac_ff, double frac_auger, double frac_tot);
double planck_spectral_radiance(double nu, void *params);
double nu_times_radiance(double nu, void *params);
double mean_frequency_nu_range(double T, double nu_min, double nu_max);
double estimate_temperature_from_mean_frequency(double mean_nu_target, double nu_min, double nu_max, double initial_guess);
int normalise_simple_estimators(PlasmaPtr xplasma);
void update_persistent_directional_flux_estimators(int nplasma, double flux_persist_scale);
/* extract.c */
int extract(WindPtr w, PhotPtr p, int itype);
int extract_one(WindPtr w, PhotPtr pp, int nspec);
/* frame.c */
int check_frame(PhotPtr p, enum frame desired_frame, char *msg);
double calculate_gamma_factor(double vel[3]);
int observer_to_local_frame(PhotPtr p_in, PhotPtr p_out);
int local_to_observer_frame(PhotPtr p_in, PhotPtr p_out);
int observer_to_local_frame_disk(PhotPtr p_in, PhotPtr p_out);
int local_to_observer_frame_disk(PhotPtr p_in, PhotPtr p_out);
double observer_to_local_frame_ds(PhotPtr p_obs, double ds_obs);
double local_to_observer_frame_ds(PhotPtr p_obs, double ds_cmf);
double observer_to_local_frame_velocity(double *v_obs, double *v, double *v_cmf);
double local_to_observer_frame_velocity(double *v_cmf, double *v, double *v_obs);
int local_to_observer_frame_ruler_transform(double v[], double dx_cmf[], double dx_obs[]);
int observer_to_local_frame_ruler_transform(double v[], double dx_obs[], double dx_cmf[]);
int lorentz_transform(PhotPtr p_in, PhotPtr p_out, double v[]);
/* gradv.c */
double dvwind_ds_cmf(PhotPtr p);
int calculate_cell_dvds_ave(int ndom, WindPtr cell);
int calculate_cell_dvds_max(int ndom, WindPtr cell);
double get_dvds_max(PhotPtr p);
/* gridwind.c */
int create_wind_and_plasma_cell_maps(void);
int calloc_wind(int nelem);
int calloc_plasma(int nelem);
int check_plasma(PlasmaPtr xplasma, char message[]);
int calloc_macro(int nelem);
int calloc_estimators(int nelem);
int calloc_dyn_plasma(int nelem);
int calloc_matom_matrix(int nelem);
void allocate_macro_matrix(double ***matrix_addr, int matrix_size);
/* homologous.c */
int get_homologous_params(int ndom);
double homologous_velocity(int ndom, double x[], double v[]);
double homologous_rho(int ndom, double x[]);
/* hydro_import.c */
int get_hydro_wind_params(int ndom);
int get_hydro(int ndom);
double hydro_velocity(int ndom, double x[], double v[]);
double hydro_rho(double x[]);
double hydro_temp(double x[]);
int rtheta_make_hydro_grid(int ndom, WindPtr w);
int rtheta_hydro_cell_volume(WindPtr w);
int hydro_frac(double coord, double coord_array[], int imax, int *cell1, int *cell2, double *frac);
double hydro_interp_value(double array[], int im, int ii, int jm, int jj, double f1, double f2);
int hydro_restart(int ndom);
void create_hydro_output_files(void);
/* import.c */
int import_wind(int ndom);
int import_wind2(int ndom, char *filename);
int import_set_wind_boundaries(int ndom);
int import_make_grid(int ndom, WindPtr w);
double import_velocity(int ndom, double *x, double *v);
double import_rho(int ndom, double *x);
double import_temperature(int ndom, double *x, int return_t_e);
/* import_calloc.c */
void calloc_import(int coord_type, int ndom);
void free_import(int coord_type, int ndom);
/* import_cylindrical.c */
int import_cylindrical(int ndom, char *filename);
int import_cylindrical_setup_boundaries(int ndom);
int cylindrical_make_grid_import(WindPtr w, int ndom);
double velocity_cylindrical(int ndom, double *x, double *v);
double rho_cylindrical(int ndom, double *x);
double temperature_cylindrical(int ndom, double *x, int return_t_e);
/* import_rtheta.c */
int import_rtheta(int ndom, char *filename);
int import_rtheta_setup_boundaries(int ndom);
int rtheta_make_grid_import(WindPtr w, int ndom);
double velocity_rtheta(int ndom, double *x, double *v);
double rho_rtheta(int ndom, double *x);
double temperature_rtheta(int ndom, double *x, int return_t_e);
/* import_spherical.c */
int import_1d(int ndom, char *filename);
int import_spherical_setup_boundaries(int ndom);
int spherical_make_grid_import(WindPtr w, int ndom);
double velocity_1d(int ndom, double *x, double *v);
double rho_1d(int ndom, double *x);
double temperature_1d(int ndom, double *x, int return_t_e);
/* ionization.c */
void update_old_plasma_variables(PlasmaPtr xplasma);
int ion_abundances(PlasmaPtr xplasma, int mode);
int convergence(PlasmaPtr xplasma);
int check_convergence(void);
int one_shot(PlasmaPtr xplasma, int mode);
double calc_te(PlasmaPtr xplasma, double tmin, double tmax);
double zero_emit(double t);
double zero_emit2(double t, void *params);
/* janitor.c */
void free_domains(void);
void free_wind_grid(void);
void free_plasma_grid(void);
void free_macro_grid(void);
void free_photons(void);
void free_atomic_data(void);
void free_spectra(void);
void clean_on_exit(void);
/* knigge.c */
int get_knigge_wind_params(int ndom);
double kn_velocity(int ndom, double x[], double v[]);
double kn_rho(int ndom, double x[]);
double kn_vzero(double r);
double kn_wind_mdot_integral(double r, void *params);
double kn_rho_zero(int ndom, double r);
/* levels.c */
int levels(PlasmaPtr xplasma, int mode);
/* lines.c */
double total_line_emission(PlasmaPtr xplasma, double f1, double f2);
double lum_lines(PlasmaPtr xplasma, int nmin, int nmax);
double two_level_atom(struct lines *line_ptr, PlasmaPtr xplasma, double *d1, double *d2);
double line_nsigma(struct lines *line_ptr, PlasmaPtr xplasma);
double scattering_fraction(struct lines *line_ptr, PlasmaPtr xplasma);
double p_escape(struct lines *line_ptr, PlasmaPtr xplasma);
double p_escape_from_tau(double tau);
int line_heat(PlasmaPtr xplasma, PhotPtr pp, int nres);
/* macro_accelerate.c */
void calc_matom_matrix(PlasmaPtr xplasma, double **matom_matrix);
int fill_kpkt_rates(PlasmaPtr xplasma, int *escape, PhotPtr p);
double f_matom_emit_accelerate(PlasmaPtr xplasma, int upper, double freq_min, double freq_max);
double f_kpkt_emit_accelerate(PlasmaPtr xplasma, double freq_min, double freq_max);
int matom_deactivation_from_matrix(PlasmaPtr xplasma, int uplvl);
int calc_all_matom_matrices(void);
/* macro_gen_f.c */
double get_matom_f(int mode);
double get_matom_f_accelerate(int mode);
/* macro_gov.c */
int macro_gov(PhotPtr p, int *nres, int matom_or_kpkt, int *which_out);
int macro_pops(PlasmaPtr xplasma, double xne);
int macro_pops_fill_rate_matrix(MacroPtr mplasma, PlasmaPtr xplasma, double xne, int index_element, double rate_matrix[600][600], int radiative_flag[600][600], int conf_to_matrix[600]);
int macro_pops_check_for_population_inversion(int index_element, double *populations, int radiative_flag[600][600], int conf_to_matrix[600]);
int macro_pops_check_densities_for_numerical_errors(PlasmaPtr xplasma, int index_element, double *populations, int conf_to_matrix[600], int n_iterations);
void macro_pops_copy_to_xplasma(PlasmaPtr xplasma, int index_element, double *populations, int conf_to_matrix[600]);
/* matom.c */
int matom(PhotPtr p, int *nres, int *escape);
double b12(struct lines *line_ptr);
double xalpha_sp(struct topbase_phot *cont_ptr, PlasmaPtr xplasma, int ichoice);
double alpha_sp(struct topbase_phot *cont_ptr, PlasmaPtr xplasma, int ichoice);
double scaled_alpha_sp_integral_band_limited(struct topbase_phot *cont_ptr, PlasmaPtr xplasma, int ichoice, double freq_min, double freq_max);
double alpha_sp_integrand(double freq, void *params);
int kpkt(PhotPtr p, int *nres, int *escape, int mode);
int fake_matom_bb(PhotPtr p, int *nres, int *escape);
int fake_matom_bf(PhotPtr p, int *nres, int *escape);
int emit_matom(WindPtr w, PhotPtr p, int *nres, int upper, double freq_min, double freq_max);
/* matom_diag.c */
int matom_emiss_report(void);
/* matrix_cpu.c */
const char *get_matrix_error_string(int error_code);
int solve_matrix(double *a_matrix, double *b_matrix, int size, double *x_matrix, int nplasma);
int invert_matrix(double *matrix, double *inverted_matrix, int num_rows);
/* matrix_ion.c */
int matrix_ion_populations(PlasmaPtr xplasma, int mode);
int populate_ion_rate_matrix(double rate_matrix[nions][nions], double pi_rates[nions], double inner_rates[n_inner_tot], double rr_rates[nions], double b_temp[nions], double xne, double nh1, double nh2);
/* matrix_ion2.c */
int matrix_ion_populations2(PlasmaPtr xplasma, int mode);
/* models_extern_init.c */
/* para_update.c */
int get_parallel_nrange(int rank, int ntotal, int nproc, int *my_nmin, int *my_nmax);
int get_max_cells_per_rank(const int n_total);
int calculate_comm_buffer_size(const int num_ints, const int num_doubles);
/* parse.c */
int parse_command_line(int argc, char *argv[]);
void help(void);
/* partition.c */
int partition_functions(PlasmaPtr xplasma, int mode);
int partition_functions_2(PlasmaPtr xplasma, int xnion, double temp, double weight);
/* paths.c */
Wind_Paths_Ptr wind_paths_constructor(WindPtr wind);
int reverb_init(WindPtr wind);
int wind_paths_init(WindPtr wind);
int line_paths_add_phot(WindPtr wind, PhotPtr pp, int *nres);
int wind_paths_add_phot(WindPtr wind, PhotPtr pp);
int simple_paths_gen_phot(PhotPtr pp);
double r_draw_from_path_histogram(Wind_Paths_Ptr PathPtr);
int wind_paths_gen_phot(WindPtr wind, PhotPtr pp);
int line_paths_gen_phot(WindPtr wind, PhotPtr pp, int nres);
int wind_paths_evaluate_single(Wind_Paths_Ptr paths);
int wind_paths_evaluate(WindPtr wind, int i_rank);
int wind_paths_dump(WindPtr wind, int rank_global);
int wind_paths_output_dump(WindPtr wind, int i_rank);
int wind_paths_point_index(int i, int j, int k, int i_top, DomainPtr dom);
int wind_paths_sphere_point_index(int i, int j, int k);
int wind_paths_output_vtk(WindPtr wind, int ndom);
/* phot_util.c */
int init_dummy_phot(PhotPtr p);
int stuff_phot(PhotPtr pin, PhotPtr pout);
int move_phot(PhotPtr pp, double ds);
int comp_phot(PhotPtr p1, PhotPtr p2);
double ds_to_cone(ConePtr cc, struct photon *p);
double ds_to_sphere(double r, struct photon *p);
double ds_to_sphere2(double x[], double r, struct photon *p);
int quadratic(double a, double b, double c, double r[]);
double ds_to_plane(struct plane *pl, struct photon *p, int force_positive_z);
double ds_to_closest_approach(double x[], struct photon *p, double *impact_parameter);
double ds_to_cylinder(double rho, struct photon *p);
/* photon2d.c */
int translate(WindPtr w, PhotPtr pp, double tau_scat, double *tau, int *nres);
int translate_in_space(PhotPtr pp);
double ds_to_wind(PhotPtr pp, int *ndom_current);
int translate_in_wind(WindPtr w, PhotPtr p, double tau_scat, double *tau, int *nres);
double smax_in_cell(PhotPtr p);
double ds_in_cell(int ndom, PhotPtr p);
/* photon_gen.c */
int define_phot(PhotPtr p, double f1, double f2, long nphot_tot, int ioniz_or_extract, int iwind, int freq_sampling);
double populate_bands(int ioniz_or_extract, int iwind, struct xbands *band);
int xdefine_phot(double f1, double f2, int ioniz_or_extract, int iwind, int print_mode, int tot_flag);
int phot_status(void);
int xmake_phot(PhotPtr p, double f1, double f2, int ioniz_or_extract, int iwind, double weight, int iphot_start, int nphotons);
int star_init(double freqmin, double freqmax, int ioniz_or_extract, double *f);
int photo_gen_star(PhotPtr p, double r, double t, double weight, double f1, double f2, int spectype, int istart, int nphot);
double bl_init(double lum_bl, double t_bl, double freqmin, double freqmax, int ioniz_or_extract, double *f);
int photon_checks(PhotPtr p, double freqmin, double freqmax, char *comment);
/* photon_gen_matom.c */
double get_kpkt_f(void);
double get_kpkt_heating_f(void);
int photo_gen_kpkt(PhotPtr p, double weight, int photstart, int nphot);
int photo_gen_matom(PhotPtr p, double weight, int photstart, int nphot);
/* pi_rates.c */
double calc_pi_rate(int nion, PlasmaPtr xplasma, int mode, int type);
double tb_planck(double freq, void *params);
double tb_logpow(double freq, void *params);
double tb_exp(double freq, void *params);
/* sirocco_extern_init.c */
/* radiation.c */
double radiation(PhotPtr p, double ds);
double kappa_ff(PlasmaPtr xplasma, double freq);
double sigma_phot(struct topbase_phot *x_ptr, double freq);
double den_config(PlasmaPtr xplasma, int nconf);
double pop_kappa_ff_array(void);
double mean_intensity(PlasmaPtr xplasma, double freq, int mode);
double mean_intensity_from_models(PlasmaPtr xplasma, double freq, int mode);
double mean_intensity_bb_estimate(double freq, double t_r, double w);
/* random.c */
int randvec(double a[], double r);
int randvcos(double lmn[], double north[]);
double vcos(double x, void *params);
int randvdipole(double lmn[], double north[]);
double vdipole(double cos_theta, void *params);
int init_rand(int seed);
void init_rng_directory(char *root, int rank);
void save_gsl_rng_state(void);
void reload_gsl_rng_state(void);
double random_number(double min, double max);
/* rdpar.c */
int opar(char filename[]);
int add_par(char filename[]);
int cpar(char filename[]);
int rdpar_init(void);
int string_process(char question[], char dummy[]);
int string_process_from_command_line(char question[], char dummy[]);
int string_process_from_file(char question[], char dummy[]);
int rdpar_store_record(char *name, char *value);
int rdpar_save(FILE *file_ptr);
int rdpar_comment(char *format, ...);
int message(char string[]);
int rdstr(char question[], char answer[]);
int rdchar(char question[], char *answer);
int rdint(char question[], int *answer);
int rdflo(char question[], float *answer);
int rddoub(char question[], double *answer);
int rdline(char question[], char answer[]);
int string2int(char *word, char *string_choices, char *string_values, char *string_answer);
int rdchoice(char question[], char answers[], char *answer);
int get_root(char root[], char total[]);
int rdpar_set_mpi_rank(int rank);
int rdpar_set_verbose(int vlevel);
int rdpar_check(void);
/* rdpar_init.c */
int init_choices(void);
int get_choices(char *question, char *choices, struct rdpar_choices *qstruct);
/* recipes.c */
double num_int(double (*func)(double, void *), double a, double b, double eps);
double zero_find(double (*func)(double, void *), double x1, double x2, double tol, int *ierr);
double find_function_minimum(double a, double m, double b, double (*func)(double, void *), double tol, double *xmin);
int fraction(double value, double array[], int npts, int *ival, double *f, int mode);
int linterp(double x, double xarray[], double yarray[], int xdim, double *y, int mode);
/* recomb.c */
double fb_topbase_partial(double freq);
double fb_topbase_partial2(double freq, void *params);
double integ_fb(double t, double f1, double f2, int nion, int fb_choice, int mode);
double total_fb(PlasmaPtr xplasma, double t, double f1, double f2, int fb_choice, int mode);
double one_fb(PlasmaPtr xplasma, double f1, double f2);
int num_recomb(PlasmaPtr xplasma, double t_e, int mode);
double fb(PlasmaPtr xplasma, double t, double freq, int ion_choice, int fb_choice);
int init_freebound(double t1, double t2, double f1, double f2);
double get_nrecomb(double t, int nion, int mode);
double get_fb(double t, int nion, int narray, int fb_choice, int mode);
double xinteg_fb(double t, double f1, double f2, int nion, int fb_choice);
double xinteg_inner_fb(double t, double f1, double f2, int nion, int fb_choice);
double total_rrate(int nion, double T);
double gs_rrate(int nion, double T);
int sort_and_compress(double *array_in, double *array_out, int npts);
int compare_doubles(const void *a, const void *b);
double matom_select_bf_freq(WindPtr one, int nconf);
/* resonate.c */
double calculate_ds(WindPtr w, PhotPtr p, double tau_scat, double *tau, int *nres, double smax, int *istat);
int select_continuum_scattering_process(double kap_cont, double kap_es, double kap_ff, PlasmaPtr xplasma);
double kappa_bf(PlasmaPtr xplasma, double freq, int macro_all);
int kbf_need(double freq_min, double freq_max);
double sobolev(WindPtr one, double x[], double den_ion, struct lines *lptr, double dvds);
int scatter(PhotPtr p, int *nres, int *nnscat);
/* reverb.c */
double delay_to_observer(PhotPtr pp);
int delay_dump_prep(int restart_stat);
int delay_dump_finish(void);
int delay_dump_combine(int i_ranks);
int delay_dump(PhotPtr p, int np);
int delay_dump_single(PhotPtr pp, int i_spec);
/* roche.c */
int binary_basics(void);
int hit_secondary(PhotPtr p);
double pillbox(PhotPtr p, double *smin, double *smax);
double phi(double s, void *params);
double dphi_ds(double s, void *params);
double roche_width(double x, void *params);
double roche2_half_width(void);
/* rtheta.c */
double rtheta_ds_in_cell(int ndom, PhotPtr p);
int rtheta_make_grid(int ndom, WindPtr w);
int rtheta_make_cones(int ndom, WindPtr w);
int rtheta_wind_complete(int ndom, WindPtr w);
int rtheta_cell_volume(WindPtr w);
int rtheta_where_in_grid(int ndom, double x[]);
int rtheta_get_random_location(int n, double x[]);
int rtheta_extend_density(int ndom, WindPtr w);
int rtheta_is_cell_in_wind(int n);
/* run.c */
int calculate_ionization(int restart_stat);
int make_spectra(int restart_stat);
/* saha.c */
int nebular_concentrations(PlasmaPtr xplasma, int mode);
int concentrations(PlasmaPtr xplasma, int mode);
int saha(PlasmaPtr xplasma, double ne, double t);
int lucy(PlasmaPtr xplasma);
int lucy_mazzali1(double nh, double t_r, double t_e, double www, int nelem, double ne, double density[], double xne, double newden[]);
int fix_concentrations(PlasmaPtr xplasma, int mode);
double get_ne(double density[]);
/* setup.c */
int init_geo(void);
int get_spectype(int yesno, char *question, int *spectype);
int init_advanced_modes(void);
int init_observers(void);
PhotPtr init_photons(void);
int init_ionization(void);
double setup_dfudge(void);
void setup_atomic_data(const char *atomic_filename);
/* setup_disk.c */
double get_disk_params(void);
/* setup_domains.c */
int get_domain_params(int ndom);
void allocate_domain_wind_coords(int ndom);
int get_wind_params(int ndom);
int setup_windcone(void);
int init_windcone(double r, double z, double dzdr, int allow_negative_dzdr, ConePtr one_windcone);
/* setup_files.c */
int init_log_and_windsave(int restart_stat);
int setup_created_files(void);
/* setup_line_transfer.c */
int get_line_transfer_mode(void);
int line_transfer_help_message(void);
/* setup_reverb.c */
int get_meta_params(void);
/* setup_star_bh.c */
double get_stellar_params(void);
int get_bl_and_agn_params(double lstar);
/* shell_wind.c */
int get_shell_wind_params(int ndom);
int shell_make_grid(int ndom, WindPtr w);
/* signal.c */
int xsignal(char *root, char *format, ...);
int d_xsignal(char *root, char *format, ...);
int xsignal_rm(char *root);
int set_max_time(char *root, double t);
int check_time(char *root);
/* spectra.c */
void spectrum_allocate(int nspec);
int spectrum_init(double f1, double f2, int nangle, double angle[], double phase[], int scat_select[], int top_bot_select[], int select_extract, double rho_select[], double z_select[], double az_select[], double r_select[]);
int spectrum_create(PhotPtr p, int nangle, int select_extract);
int spec_add_one(PhotPtr p, int spec_type);
int spectrum_summary(char filename[], int nspecmin, int nspecmax, int select_spectype, double renorm, int loglin, int iwind);
int spectrum_restart_renormalise(int nangle);
/* spectral_estimators.c */
int spectral_estimators(PlasmaPtr xplasma);
double pl_alpha_func_log(double alpha);
double pl_alpha_func_log2(double alpha, void *params);
double pl_logmean(double alpha, double lnumin, double lnumax);
double pl_log_w(double j, double alpha, double lnumin, double lnumax);
double pl_log_stddev(double alpha, double lnumin, double lnumax);
double exp_temp_func(double exp_temp);
double exp_temp_func2(double exp_temp, void *params);
double exp_mean(double exp_temp, double numin, double numax);
double exp_w(double j, double exp_temp, double numin, double numax);
double exp_stddev(double exp_temp, double numin, double numax);
/* spherical.c */
double spherical_ds_in_cell(int ndom, PhotPtr p);
int spherical_make_grid(int ndom, WindPtr w);
int spherical_wind_complete(int ndom, WindPtr w);
int spherical_cell_volume(WindPtr w);
int spherical_where_in_grid(int ndom, double x[]);
int spherical_get_random_location(int n, double x[]);
int spherical_extend_density(int ndom, WindPtr w);
/* stellar_wind.c */
int get_stellar_wind_params(int ndom);
double stellar_velocity(int ndom, double x[], double v[]);
double stellar_rho(int ndom, double x[]);
/* sv.c */
int get_sv_wind_params(int ndom);
double sv_velocity(double x[], double v[], int ndom);
double sv_rho(int ndom, double x[]);
double sv_find_wind_rzero(int ndom, double p[]);
int sv_zero_init(double p[]);
double sv_zero_r(double r, void *params);
double sv_theta_wind(int ndom, double r);
double sv_wind_mdot_integral(double r, void *params);
/* synonyms.c */
int get_question_name_length(char question[]);
int are_synonym_lists_valid(void);
int is_input_line_synonym_for_question(char question[], char input_line[]);
/* time.c */
double timer(void);
int get_time(char curtime[]);
struct timeval init_timer_t0(void);
void print_timer_duration(char *msg, struct timeval timer_t0);
/* trans_phot.c */
int trans_phot(WindPtr w, PhotPtr p, int iextract);
int trans_phot_single(WindPtr w, PhotPtr p, int iextract);
/* vvector.c */
double dot(double a[], double b[]);
double length(double a[]);
int renorm(double a[], double scalar);
int rescale(double a[], double scalar, double b[]);
int cross(double a[], double b[], double c[]);
int vmove(double u[], double lmn[], double s, double result[]);
int vsub(double u[], double v[], double result[]);
int vadd(double u[], double v[], double result[]);
int stuff_v(double vin[], double vout[]);
double dot_tensor_vec(double tensor[3][3], double vin[3], double vout[3]);
int project_from_xyz_cyl(double a[], double b[], double result[]);
int project_from_cyl_xyz(double a[], double b[], double result[]);
int create_basis(double u[], double v[], struct basis *basis_new);
int project_from(struct basis *basis_from, double v_in[], double v_out[]);
int project_to(struct basis *basis_from, double v_in[], double v_out[]);
int reorient(struct basis *basis_from, struct basis *basis_to, double v_from[], double v_to[]);
/* walls.c */
int walls(PhotPtr p, PhotPtr pold, double *normal);
/* wind.c */
int where_in_wind(double x[], int *ndomain);
double model_velocity(int ndom, double x[], double v[]);
int model_vgrad(int ndom, double x[], double v_grad[][3]);
double get_div_v_in_cmf_frame(int ndom, double *x);
double model_rho(int ndom, double x[]);
int wind_check(WindPtr www, int n);
/* wind2d.c */
int where_in_grid(int ndom, double x[]);
int vwind_xyz(int ndom, PhotPtr p, double v[]);
int wind_div_v(int ndom, WindPtr cell);
double rho(WindPtr w, double x[]);
int mdot_wind(WindPtr w, double z, double rmax);
int get_random_location(int n, double x[]);
int zero_scatters(void);
int check_corners_inwind(int n);
int check_grid(void);
/* wind_sum.c */
int xtemp_rad(WindPtr w);
/* wind_updates2d.c */
int wind_update(WindPtr w);
int report_bf_simple_ionpool(void);
void wind_rad_init(void);
void check_heating_rates_for_plasma_cell(const int n_plasma);
void init_plasma_rad_properties(void);
void init_macro_rad_properties(void);
void shell_output_wind_update_diagnostics(double xsum, double psum, double fsum, double csum, double icsum, double lsum, double ausum, double chexsum, double cool_sum, double lum_sum);
/* wind_util.c */
int coord_fraction(int ndom, int ichoice, double x[], int ii[], double frac[], int *nelem);
int where_in_2dcell(int ichoice, double x[], int n, double *fx, double *fz);
int wind_n_to_ij(int ndom, int n, int *i, int *j);
int wind_ij_to_n(int ndom, int i, int j, int *n);
int wind_x_to_n(double x[], int *n);
/* windsave.c */
int wind_save(char filename[]);
int wind_read(char filename[]);
void wind_complete(void);
int spec_save(char filename[]);
int spec_read(char filename[]);
/* windsave2table_sub.c */
int do_windsave2table(char *root, int ion_switch, int edge_switch);
int create_master_table(int ndom, char rootname[]);
int create_heat_table(int ndom, char rootname[]);
int create_convergence_table(int ndom, char rootname[]);
int create_velocity_gradient_table(int ndom, char rootname[]);
int create_ion_table(int ndom, char rootname[], int iz, int ion_switch);
double *get_ion(int ndom, int element, int istate, int iswitch, char *name);
double *get_one(int ndom, char variable_name[]);
int get_one_array_element(int ndom, char variable_name[], int array_dim, double xval[]);
int create_spec_table(int ndom, char rootname[]);
int create_detailed_cell_spec_table(int ncell, char rootname[]);
int create_big_detailed_spec_table(int ndom, char *rootname);
/* xlog.c */
int Log_init(char *filename);
int Log_append(char *filename);
void Log_close(void);
int Log_set_verbosity(int vlevel);
int Log_print_max(int print_max);
int Log_quit_after_n_errors(int n);
int Log(char *format, ...);
int Log_silent(char *format, ...);
int Error(char *format, ...);
int Error_silent(char *format, ...);
int Shout(char *format, ...);
int sane_check(double x);
int error_count(char *format);
int error_summary(char *message);
int error_summary_parallel(char *msg);
int Log_flush(void);
int Log_set_mpi_rank(int rank, int n_mpi);
int Log_parallel(char *format, ...);
int Debug(char *format, ...);
void Exit(int error_code);
/* xtest.c */
int xtest(void);
/* zeta.c */
double compute_zeta(double temp, int nion, int mode);
/* swind_sub.c */
int zoom(int direction);
int overview(WindPtr w, char rootname[]);
int position_summary(WindPtr w);
int abs_summary(WindPtr w, char rootname[], int ochoice);
int shock_heating_summary(WindPtr w, char rootname[], int ochoice);
int adiabatic_cooling_summary(WindPtr w, char rootname[], int ochoice);
int lum_summary(WindPtr w, char rootname[], int ochoice);
int photo_summary(WindPtr w, char rootname[], int ochoice);
int recomb_summary(WindPtr w, char rootname[], int ochoice);
int electron_summary(WindPtr w, char rootname[], int ochoice);
int rho_summary(WindPtr w, char rootname[], int ochoice);
int plasma_cell(WindPtr w, char rootname[], int ochoice);
int freq_summary(WindPtr w, char rootname[], int ochoice);
int nphot_summary(WindPtr w, char rootname[], int ochoice);
int temp_summary(WindPtr w, char rootname[], int ochoice);
int temp_rad(WindPtr w, char rootname[], int ochoice);
int weight_summary(WindPtr w, char rootname[], int ochoice);
int velocity_summary(WindPtr w, char rootname[], int ochoice);
int mo_summary(WindPtr w, char rootname[], int ochoice);
int vol_summary(WindPtr w, char rootname[], int ochoice);
int wind_element(WindPtr w);
int tau_h_summary(WindPtr w, char rootname[], int ochoice);
int coolheat_summary(WindPtr w, char rootname[], int ochoice);
int complete_file_summary(WindPtr w, char root[], int ochoice);
int wind_reg_summary(WindPtr w, char rootname[], int ochoice);
int dvds_summary(WindPtr w, char rootname[], int ochoice);
int IP_summary(WindPtr w, char rootname[], int ochoice);
int alpha_summary(WindPtr w, char rootname[], int ochoice);
int J_summary(WindPtr w, char rootname[], int ochoice);
int J_scat_summary(WindPtr w, char rootname[], int ochoice);
int phot_split(WindPtr w, char rootname[], int ochoice);
int thompson(WindPtr w, char rootname[], int ochoice);
int nscat_split(WindPtr w, char rootname[], int ochoice);
int convergence_summary(WindPtr w, char rootname[], int ochoice);
int convergence_all(WindPtr w, char rootname[], int ochoice);
int model_bands(WindPtr w, char rootname[], int ochoice);
int heatcool_summary(WindPtr w, char rootname[], int ochoice);
int complete_physical_summary(WindPtr w, char rootname[], int ochoice);
int complete_ion_summary(WindPtr w, char rootname[], int ochoice);
double get_density_or_frac(PlasmaPtr xplasma, int element, int istate, int frac_choice);
int find_ion(int element, int istate);
int find_element(int element);
int get_los_dvds(WindPtr w, char rootname[], int ochoice);
int grid_summary(WindPtr w, char rootname[], int ochoice);
int flux_summary(WindPtr w, char rootname[], int ochoice);
/* swind_ion.c */
int ion_summary(WindPtr w, int element, int istate, int iswitch, char rootname[], int ochoice);
int tau_ave_summary(WindPtr w, int element, int istate, double freq, char rootname[], int ochoice);
int line_summary(WindPtr w, char rootname[], int ochoice);
int total_emission_summary(char rootname[], int ochoice);
int modify_te(WindPtr w, char rootname[], int ochoice);
int partial_measure_summary(WindPtr w, int element, int istate, char rootname[], int ochoice);
int collision_summary(WindPtr w, char rootname[], int ochoice);
/* swind_write.c */
int write_array(char filename[], int choice);
int display(char name[]);
/* swind_macro.c */
int xadiabatic_cooling_summary(WindPtr w, char rootname[], int ochoice);
int macro_summary(WindPtr w, char rootname[], int ochoice);
int ion_overview(int icell);
int config_overview(int n, int icell);
int depcoef_overview(int icell);
int copy_plasma(PlasmaPtr x1, PlasmaPtr x2);
int dealloc_copied_plasma(PlasmaPtr xcopy);
int depcoef_overview_specific(int version, int nconfig, WindPtr w, char rootname[], int ochoice);
int level_popsoverview(int nplasma, WindPtr w, char rootname[], int ochoice);
int level_emissoverview(int nlev, WindPtr w, char rootname[], int ochoice);
int level_escapeoverview(int nlev, WindPtr w, char rootname[], int ochoice);
int level_tauoverview(int nlev, WindPtr w, char rootname[], int ochoice);
/* swind.c */
int main(int argc, char *argv[]);
int one_choice(int choice, char *root, int ochoice);
void swind_help(void);
/* windsave2table.c */
void parse_arguments(int argc, char *argv[], char root[], int *ion_switch, int *spec_switch, int *edge_switch);
int main(int argc, char *argv[]);
/* windsave2table_sub.c */
int do_windsave2table(char *root, int ion_switch, int edge_switch);
int create_master_table(int ndom, char rootname[]);
int create_heat_table(int ndom, char rootname[]);
int create_convergence_table(int ndom, char rootname[]);
int create_velocity_gradient_table(int ndom, char rootname[]);
int create_ion_table(int ndom, char rootname[], int iz, int ion_switch);
double *get_ion(int ndom, int element, int istate, int iswitch, char *name);
double *get_one(int ndom, char variable_name[]);
int get_one_array_element(int ndom, char variable_name[], int array_dim, double xval[]);
int create_spec_table(int ndom, char rootname[]);
int create_detailed_cell_spec_table(int ncell, char rootname[]);
int create_big_detailed_spec_table(int ndom, char *rootname);
