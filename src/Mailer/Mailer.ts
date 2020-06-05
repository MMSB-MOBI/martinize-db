import nodemailer from 'nodemailer';
import { TEMPLATE_DIR, URLS, DEFAULT_MAILER_NAME, DEFAULT_MAILER_ADDRESS, MAILER_ENFORCE_RECIPIENT, MAILER_TRANSPORT_SETTINGS } from '../constants';
import Twig from 'twig';
import logger from '../logger';

export default new class Mailer {
  protected transporter = nodemailer.createTransport(MAILER_TRANSPORT_SETTINGS);

  public default_sender = { name: DEFAULT_MAILER_NAME, address: DEFAULT_MAILER_ADDRESS };

  async send(send_options: nodemailer.SendMailOptions, template_name: string, options: { [variableName: string]: any }) {
    if (!send_options.to) {
      throw new Error("You must define a mail recipient.");
    }

    if (!options.site_url) {
      options.site_url = URLS.SERVER;
    }

    const file = TEMPLATE_DIR + template_name + (template_name.endsWith('.twig') ? "" : ".twig");

    const content = await new Promise((resolve, reject) => {
      // @ts-ignore Incorrect typedef for options
      Twig.renderFile(file, options, (err: Error, res: string) => {
        if (err) {
          reject(err);
        }
        resolve(res);
      })
    }) as string;

    send_options.html = content;

    if (!send_options.from) {
      send_options.from = this.default_sender;
      send_options.sender = this.default_sender;
    }
    if (!send_options.subject && options.title) {
      send_options.subject = options.title;
    }

    if (MAILER_ENFORCE_RECIPIENT) {
      send_options.to = MAILER_ENFORCE_RECIPIENT;
    }

    return this.mail(send_options);
  }

  protected async mail(options: nodemailer.SendMailOptions) {
    const info = await this.transporter.sendMail(options);

    return info as { messageId: string };
  }
}();
